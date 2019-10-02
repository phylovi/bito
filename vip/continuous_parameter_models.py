import numpy as np
import pandas as pd
import scipy as sp

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp

tf.enable_v2_behavior()
tfd = tfp.distributions


def exponential_factory(params):
    return tfd.Exponential(rate=params[:, 0])


def gamma_factory(params):
    return tfd.Gamma(concentration=tf.exp(params[:, 0]), rate=tf.exp(params[:, 1]))


def lognormal_factory(params):
    return tfd.LogNormal(loc=params[:, 0], scale=params[:, 1])


def truncated_lognormal_factory(params):
    exp_shift = tfp.bijectors.Chain(
        [
            tfp.bijectors.AffineScalar(shift=-tf.math.exp(params[:, 2])),
            tfp.bijectors.Exp(),
        ]
    )
    return tfd.TransformedDistribution(
        distribution=tfd.TruncatedNormal(
            loc=params[:, 0], scale=params[:, 1], low=params[:, 2], high=999
        ),
        bijector=exp_shift,
        name="TruncatedLogNormal",
    )


class TFContinuousParameterModel:
    """An object to model a collection of variables with a given q
    distribution.

    Each of these variables uses a number of parameters which we
    optimize using SGD variants.

    Samples are laid out as particles x variables, which we will call z-shape.

    * log_p is the target probability.
    """

    def __init__(self, q_factory, initial_params, variable_count, particle_count):
        assert initial_params.ndim == 1
        self.q_factory = q_factory
        self.name = self.q_factory(np.array([initial_params])).name
        self.param_matrix = np.full(
            (variable_count, len(initial_params)), initial_params
        )
        self.particle_count = particle_count
        # The current stored sample.
        self.z = None
        # The gradient of z with respect to the parameters of q.
        self.grad_z = None
        # The stochastic gradient of log sum q for z.
        self.grad_log_sum_q = None

    @property
    def variable_count(self):
        return self.param_matrix.shape[0]

    @property
    def param_count(self):
        return self.param_matrix.shape[1]

    def mode_match(self, modes):
        """Some crazy heuristics for mode matching with the given branch
        lengths."""
        log_modes = np.log(np.clip(modes, 1e-6, None))
        biclipped_log_modes = np.log(np.clip(modes, 1e-6, 1-1e-6))
        # TODO do we want to have the lognormal variance be in log space? Or do we want
        # to clip it?
        if self.name == "LogNormal":
            self.param_matrix[:, 1] = -0.1 * biclipped_log_modes
            self.param_matrix[:, 0] = np.square(self.param_matrix[:, 1]) + log_modes
        elif self.name == "TruncatedLogNormal":
            self.param_matrix[:, 1] = -0.1 * biclipped_log_modes
            self.param_matrix[:, 0] = np.square(self.param_matrix[:, 1]) + log_modes
            self.param_matrix[:, 2] = -5
        elif self.name == "Gamma":
            self.param_matrix[:, 1] = np.log(-60.0 * biclipped_log_modes)
            self.param_matrix[:, 0] = np.log(1 + modes * self.param_matrix[:, 1])
        else:
            print("Mode matching not implemented for " + self.name)

    def suggested_step_size(self):
        return np.average(np.abs(self.param_matrix), axis=0) / 100

    def sample(self, particle_count):
        q_distribution = self.q_factory(self.param_matrix)
        return q_distribution.sample(particle_count).numpy()

    def sample_and_prep_gradients(self):
        with tf.GradientTape(persistent=True) as g:
            tf_params = tf.constant(self.param_matrix, dtype=tf.float32)
            g.watch(tf_params)
            q_distribution = self.q_factory(tf_params)
            tf_z = q_distribution.sample(self.particle_count)
            q_term = tf.math.reduce_sum(tf.math.log(q_distribution.prob(tf_z)))
        self.z = tf_z.numpy()
        # The Jacobian is laid out as particles x edges x edges x params.
        self.grad_z = np.sum(g.jacobian(tf_z, tf_params).numpy(), axis=2)
        self.grad_log_sum_q = g.gradient(q_term, tf_params).numpy()
        del g  # Should happen anyway but being explicit to remember.
        return self.z

    def clear_sample(self):
        self.z = None
        self.grad_z = None
        self.grad_log_sum_q = None

    def elbo_estimate(self, log_p, particle_count=None):
        """A naive Monte Carlo estimate of the ELBO.

        log_p must take an argument of z-shape.
        """
        if particle_count is None:
            particle_count = self.particle_count
        q_distribution = self.q_factory(self.param_matrix)
        z = q_distribution.sample(particle_count)
        return (
            np.sum(log_p(z) - np.sum(q_distribution.log_prob(z), axis=1))
        ) / particle_count

    @staticmethod
    def _chain_rule(grad_log_p_z, grad_z):
        return (
            np.tensordot(grad_log_p_z.transpose(), grad_z, axes=1)
            .diagonal()
            .transpose()
        )

    @staticmethod
    def _slow_chain_rule(grad_log_p_z, grad_z):
        particle_count, variable_count, param_count = grad_z.shape
        result = np.zeros((variable_count, param_count))
        for variable in range(variable_count):
            for param in range(param_count):
                for particle in range(particle_count):
                    result[variable, param] += (
                        grad_log_p_z[particle, variable]
                        * grad_z[particle, variable, param]
                    )
        return result

    def elbo_gradient_using_current_sample(self, grad_log_p_z):
        assert self.grad_z is not None
        if not np.allclose(
            self._chain_rule(grad_log_p_z, self.grad_z),
            self._slow_chain_rule(grad_log_p_z, self.grad_z),
        ):
            print("Warning: chain rule isn't close. Here's the difference:")
            print(
                self._chain_rule(grad_log_p_z, self.grad_z)
                - self._slow_chain_rule(grad_log_p_z, self.grad_z)
            )
        unnormalized_result = (
            self._chain_rule(grad_log_p_z, self.grad_z) - self.grad_log_sum_q
        )
        return unnormalized_result / self.particle_count

    def gradient_step(self, optimizer, step_size, grad_log_p_z, history=None):
        """Return if the gradient step was successful."""
        grad = self.elbo_gradient_using_current_sample(grad_log_p_z)
        if not np.isfinite(np.array([grad])).all():
            self.clear_sample()
            return False
        update_dict = optimizer.adam(
            {"params": step_size}, {"params": self.param_matrix}, {"params": grad}
        )
        self.param_matrix += update_dict["params"]
        self.clear_sample()
        if history is not None:
            history.append(self.param_matrix.copy())
        return True

    def likelihoods_on_grid(self, num=50, max_x=0.5):
        min_x = max_x / 100
        z = np.zeros((num, self.variable_count))
        x_vals = np.linspace(min_x, max_x, num)
        for variable in range(self.variable_count):
            z[:, variable] = x_vals
        q_distribution = self.q_factory(self.param_matrix)
        df = pd.DataFrame(
            np.apply_along_axis(
                sp.special.softmax, 0, q_distribution.log_prob(z).numpy()
            )
        )
        df["x"] = x_vals
        return df.melt(id_vars="x")
