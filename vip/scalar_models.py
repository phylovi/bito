import abc
import numpy as np

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp

tf.enable_v2_behavior()
tfd = tfp.distributions


class ScalarModel(abc.ABC):
    """An abstract base class for Scalar Models.

    See the tex file in doc to understand the notation here.

    A "variable" here is some scalar variable that's part of our model, and we will use
    "theta" to denote it. For an unrooted tree, the variables are the branch lengths for
    the ensemble of trees in the SBN support.

    A full PSP parameterization (issue #119) we will need to generalize this.

    * theta_sample is the current stored sample and is laid out as (particle, variable).
    * dg_dpsi is the gradient of the reparametrization function with respect to the
      scalar model parameters, and is laid out as (particle, variable, parameter).
    * dlog_sum_q_dpsi is the gradient of the sum of the q function over particles with
      respect to the scalar model parameters, and is laid out as (variable, parameter).

    * p is the target probability.
    * q is the distribution that we're fitting to p.
    """

    def __init__(self, initial_params, variable_count, particle_count):
        assert initial_params.ndim == 1
        self.q_params = np.full((variable_count, len(initial_params)), initial_params)
        self.particle_count = particle_count
        self.theta_sample = None
        self.dg_dpsi = None
        self.dlog_sum_q_dpsi = None

    @property
    def variable_count(self):
        return self.q_params.shape[0]

    @property
    def param_count(self):
        return self.q_params.shape[1]

    def clear_sample(self):
        self.theta_sample = None
        self.dg_dpsi = None
        self.dlog_sum_q_dpsi = None

    def suggested_step_size(self):
        return np.average(np.abs(self.q_params), axis=0) / 100

    @abc.abstractmethod
    def mode_match(self, modes):
        pass

    @abc.abstractmethod
    def sample(self, particle_count, which_variables, log_prob=False):
        pass

    @abc.abstractmethod
    def sample_and_prep_gradients(self, which_variables):
        """Sample the variables in which_variables and prepare so we can take a
        gradient with respect to them."""
        pass


class LogNormalModel(ScalarModel):
    """A log-normal model with hand-computed gradients."""

    def __init__(self, initial_params, variable_count, particle_count):
        super().__init__(initial_params, variable_count, particle_count)
        self.name = "LogNormal"

    def mu(self, which_variables):
        return self.q_params[which_variables, 0]

    def sigma(self, which_variables):
        return self.q_params[which_variables, 1]

    def mode_match(self, modes):
        """Some crazy heuristics for mode matching with the given theta
        values."""
        log_modes = np.log(np.clip(modes, 1e-6, None))
        biclipped_log_modes = np.log(np.clip(modes, 1e-6, 1 - 1e-6))
        # Issue #118: do we want to have the lognormal variance be in log space? Or do
        # we want to clip it?
        self.q_params[:, 1] = -0.1 * biclipped_log_modes
        self.q_params[:, 0] = (
            np.square(self.sigma(np.arange(self.variable_count))) + log_modes
        )

    def sample(self, particle_count, which_variables):
        return np.random.lognormal(
            self.mu(which_variables),
            self.sigma(which_variables),
            (particle_count, len(which_variables)),
        )

    def log_prob(self, theta_sample, which_variables):
        """Return a log probability for each of the particles given in
        theta_sample."""
        assert theta_sample.shape[1] == which_variables.size
        log_z = np.log(theta_sample)
        # Here's the fancy stable version.
        # ratio = np.exp(np.log((np.log(np.log(theta_sample)-mu)**2) - np.log(sigma**2))
        ratio = (log_z - self.mu(which_variables)) ** 2 / (
            2 * self.sigma(which_variables) ** 2
        )
        # Below we're summing over axis 1, which is the variable axis.
        result = (
            -0.5 * np.log(2 * np.pi)
            - np.sum(log_z, axis=1)
            - np.sum(np.log(self.sigma(which_variables)))
            - np.sum(ratio, axis=1)
        )
        return result

    def sample_and_prep_gradients(self, which_variables):
        self.theta_sample = self.sample(self.particle_count, which_variables)
        epsilon = (np.log(self.theta_sample) - self.mu(which_variables)) / self.sigma(
            which_variables
        )
        # See the tex for details about this gradient.
        self.dg_dpsi = np.empty((self.particle_count, self.variable_count, 2))
        self.dg_dpsi[:, which_variables, 0] = self.theta_sample
        self.dg_dpsi[:, which_variables, 1] = self.theta_sample * epsilon
        self.dlog_sum_q_dpsi = np.empty((self.variable_count, 2))
        self.dlog_sum_q_dpsi[which_variables, 0] = -1.0 * self.particle_count
        self.dlog_sum_q_dpsi[which_variables, 1] = -np.sum(
            epsilon, axis=0
        ) - self.particle_count / self.sigma(which_variables)


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


class TFScalarModel(ScalarModel):
    """An object to model a collection of scalar variables with a given
    distribution type via TensorFlow.

    See ScalarModel for more information.
    """

    def __init__(self, q_factory, initial_params, variable_count, particle_count):
        super().__init__(initial_params, variable_count, particle_count)
        self.q_factory = q_factory
        self.name = "TF" + self.q_factory(np.array([initial_params])).name

    def mode_match(self, modes):
        """Some crazy heuristics for mode matching with the given branch
        lengths."""
        log_modes = np.log(np.clip(modes, 1e-6, None))
        biclipped_log_modes = np.log(np.clip(modes, 1e-6, 1 - 1e-6))
        # Issue #118: do we want to have the lognormal variance be in log space? Or do
        # we want to clip it?
        if self.name == "TFLogNormal":
            self.q_params[:, 1] = -0.1 * biclipped_log_modes
            self.q_params[:, 0] = np.square(self.q_params[:, 1]) + log_modes
        elif self.name == "TFTruncatedLogNormal":
            self.q_params[:, 1] = -0.1 * biclipped_log_modes
            self.q_params[:, 0] = np.square(self.q_params[:, 1]) + log_modes
            self.q_params[:, 2] = -5
        elif self.name == "TFGamma":
            self.q_params[:, 1] = np.log(-60.0 * biclipped_log_modes)
            self.q_params[:, 0] = np.log(1 + modes * self.q_params[:, 1])
        else:
            print("Mode matching not implemented for " + self.name)

    def sample(self, particle_count, log_prob=False):
        q_distribution = self.q_factory(self.q_params)
        theta_sample = q_distribution.sample(particle_count).numpy()
        if log_prob:
            return theta_sample, np.sum(q_distribution.log_prob(theta_sample), axis=1)
        else:
            return theta_sample

    def sample_and_prep_gradients(self, which_variables):
        """Take a sample from q and prepare a gradient of the sample and of log
        q with respect to the parameters of q."""
        raise NotImplementedError("We don't use which_variables yet.")
        with tf.GradientTape(persistent=True) as g:
            tf_params = tf.constant(self.q_params, dtype=tf.float32)
            g.watch(tf_params)
            q_distribution = self.q_factory(tf_params)
            tf_z = q_distribution.sample(self.particle_count)
            q_term = tf.math.reduce_sum(q_distribution.log_prob(tf_z))
        self.theta_sample = tf_z.numpy()
        # The Jacobian is laid out as particles x edges x edges x params.
        self.dg_dpsi = np.sum(g.jacobian(tf_z, tf_params).numpy(), axis=2)
        self.dlog_sum_q_dpsi = g.gradient(q_term, tf_params).numpy()
        del g  # Should happen anyway but being explicit to remember.
        return self.theta_sample


def of_name(name, *, variable_count, particle_count):
    def build_tf_model(q_factory, initial_params):
        return TFScalarModel(
            q_factory, np.array(initial_params), variable_count, particle_count
        )

    if name == "lognormal":
        return LogNormalModel(np.array([-2.0, 0.5]), variable_count, particle_count)
    # if name == "tf_lognormal":
    #     return build_tf_model(lognormal_factory, [-2.0, 0.5])
    # if name == "tf_truncated_lognormal":
    #     return build_tf_model(truncated_lognormal_factory, [-1.0, 0.5, 0.1])
    # if name == "tf_gamma":
    #     return build_tf_model(gamma_factory, [1.3, 3.0])

    raise Exception(f"Model {name} not known.")
