import abc
import numpy as np

import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp

tf.enable_v2_behavior()
tfd = tfp.distributions


class ScalarModel(abc.ABC):
    """An abstract base class for Scalar Models.

    See the tex file in `doc` to understand the notation here, and Burrito for notes
    about abbreviations.

    * p is the target probability.
    * q is the distribution that we're fitting to p.
    """

    def __init__(self, initial_params, variable_count):
        assert initial_params.ndim == 1
        self.q_params = np.full((variable_count, len(initial_params)), initial_params)

    @property
    def variable_count(self):
        return self.q_params.shape[0]

    @property
    def param_count(self):
        return self.q_params.shape[1]

    def suggested_step_size(self):
        return np.average(np.abs(self.q_params), axis=0) / 100

    @abc.abstractmethod
    def mode_match(self, modes):
        pass

    @abc.abstractmethod
    def sample(self, particle_count, which_variables):
        pass

    @abc.abstractmethod
    def sample_and_gradients(self, particle_count, px_which_variables):
        """Sample the variables in which_variables and take a gradient with
        respect to them.

        px_which_variables is a list of arrays, the ith entry of which is what
        variables we use for the ith particle, and get out the sample and some
        gradients as described above.

        On output:

        * sample is the desired sample.
        * dg_dpsi is the gradient of the reparametrization function with respect to the
          scalar model parameters, and is laid out as (particle, variable, parameter).
        * dlog_qg_dpsi is the gradient of the sum of the q function over particles
          with respect to the scalar model parameters, and is laid out as (variable,
          parameter).

        For both of the last two arrays, the variable indexing is in terms of the
        complete set of variables.

        See the tex for details.
        """
        pass

    @abc.abstractmethod
    def log_prob(self, values, which_variables):
        """Get the log probability for the given values of the variables
        specified in which_variables."""
        pass


class LogNormalModel(ScalarModel):
    """A log-normal model with hand-computed gradients."""

    def __init__(self, initial_params, variable_count):
        super().__init__(initial_params, variable_count)
        self.name = "LogNormal"

    def mu(self, which_variables=None):
        if which_variables is None:
            return self.q_params[:, 0]
        return self.q_params[which_variables, 0]

    def sigma(self, which_variables=None):
        if which_variables is None:
            return self.q_params[:, 1]
        return self.q_params[which_variables, 1]

    def mode_match(self, modes):
        """Some crazy heuristics for mode matching."""
        log_modes = np.log(np.clip(modes, 1e-6, None))
        biclipped_log_modes = np.log(np.clip(modes, 1e-6, 1 - 1e-6))
        # Issue #118: do we want to have the lognormal variance be in log space? Or do
        # we want to clip it?
        self.q_params[:, 1] = -0.1 * biclipped_log_modes
        self.q_params[:, 0] = np.square(self.sigma()) + log_modes

    def sample(self, particle_count, px_which_variables):
        """Get a sample of size particle_count from the model.

        px_which_variables: a list of arrays, the ith entry of which is what
        variables we use for the ith particle, and get out the sample. If it is None,
        sample all the variables.
        """
        if px_which_variables is None:
            return np.random.lognormal(
                self.mu(), self.sigma(), (particle_count, self.variable_count)
            )
        # else:
        # We check that px_which_variables is a fixed width in the loop below.
        which_variables_size = px_which_variables[0].size
        sample = np.empty((particle_count, which_variables_size))
        for particle_idx, which_variables in enumerate(px_which_variables):
            assert which_variables_size == which_variables.size
            sample[particle_idx, :] = np.random.lognormal(
                self.mu(which_variables), self.sigma(which_variables)
            )
        return sample

    def sample_and_gradients(self, particle_count, px_which_variables):
        # Comments of the form eq:XX refer to equations in the tex.
        # We check that px_which_variables is a fixed width in the loop below.
        which_variables_size = px_which_variables[0].size
        sample = np.empty((particle_count, which_variables_size))
        dg_dpsi = np.zeros((particle_count, self.variable_count, 2))
        dlog_qg_dpsi = np.zeros((particle_count, self.variable_count, 2))
        # eq:dlogqgdPsi
        dlog_qg_dpsi[:, :, 0] = -1.0
        for particle_idx, which_variables in enumerate(px_which_variables):
            mu = self.mu(which_variables)
            sigma = self.sigma(which_variables)
            assert which_variables_size == which_variables.size
            sample[particle_idx, :] = np.random.lognormal(mu, sigma)
            # eq:gLogNorm
            epsilon = (np.log(sample[particle_idx, :]) - mu) / sigma
            # eq:dgdPsi
            dg_dpsi[particle_idx, which_variables, 0] = sample[particle_idx, :]
            dg_dpsi[particle_idx, which_variables, 1] = (
                sample[particle_idx, :] * epsilon
            )
            # eq:dlogqgdPsi
            dlog_qg_dpsi[particle_idx, which_variables, 1] = -epsilon - 1.0 / sigma
        return (sample, dg_dpsi, dlog_qg_dpsi)

    def log_prob(self, values, which_variables):
        """Return the log joint probability of the given values.

        The density is:
        frac{1}{x sigma sqrt{2 pi}} exp(-frac{(log x - mu)^2}{2 sigma^2})

        So the log density is
        -(log x + log sigma + 0.5 log(2 pi) + frac{(log x - mu)^2}{2 sigma^2})
        """
        assert values.size == which_variables.size
        mu = self.mu(which_variables)
        sigma = self.sigma(which_variables)
        log_values = np.log(values)
        # Here's the fancy stable version.
        # ratio = np.exp(np.log((np.log(np.log(values)-mu)**2) - np.log(sigma**2))
        ratio = (log_values - mu) ** 2 / (2 * sigma ** 2)
        return -(
            np.sum(log_values)
            + np.sum(np.log(sigma))
            + which_variables.size * 0.5 * np.log(2 * np.pi)
            + np.sum(ratio)
        )


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

    def __init__(self, q_factory, initial_params, variable_count):
        super().__init__(initial_params, variable_count)
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

    def sample(self, particle_count, px_which_variables):
        if px_which_variables is None:
            q_distribution = self.q_factory(self.q_params)
            return q_distribution.sample(particle_count).numpy()
        # else:
        # We check that px_which_variables is a fixed width in the loop below.
        which_variables_size = px_which_variables[0].size
        sample = np.empty((particle_count, which_variables_size))
        for particle_idx, which_variables in enumerate(px_which_variables):
            assert which_variables_size == which_variables.size
            q_distribution = self.q_factory(self.q_params[which_variables, :])
            sample[particle_idx, :] = q_distribution.sample().numpy()
        return sample

    def sample_and_gradients(self, particle_count, px_which_variables):
        which_variables_size = px_which_variables[0].size
        sample = np.empty((particle_count, which_variables_size))
        dg_dpsi = np.zeros((particle_count, self.variable_count, 2))
        dlog_qg_dpsi = np.zeros((particle_count, self.variable_count, 2))
        for particle_idx, which_variables in enumerate(px_which_variables):
            with tf.GradientTape(persistent=True) as g:
                tf_params = tf.constant(
                    self.q_params[which_variables, :], dtype=tf.float32
                )
                g.watch(tf_params)
                q_distribution = self.q_factory(tf_params)
                tf_sample = q_distribution.sample()
                q_term = q_distribution.log_prob(tf_sample)
            sample[particle_idx, :] = tf_sample.numpy()
            dg_dpsi[particle_idx, which_variables, :] = g.gradient(
                tf_sample, tf_params
            ).numpy()
            dlog_qg_dpsi[particle_idx, which_variables, :] = g.gradient(
                q_term, tf_params
            ).numpy()
            del g  # Should happen anyway but being explicit to remember.
        return (sample, dg_dpsi, dlog_qg_dpsi)

    def log_prob(self, values, which_variables):
        """Return the log joint probability of the given values."""
        assert values.size == which_variables.size
        return np.sum(
            self.q_factory(self.q_params[which_variables, :]).log_prob(values).numpy()
        )


def of_name(name, *, variable_count):
    def build_tf_model(q_factory, initial_params):
        return TFScalarModel(q_factory, np.array(initial_params), variable_count)

    if name == "lognormal":
        return LogNormalModel(np.array([-2.0, 0.5]), variable_count)
    if name == "tf_lognormal":
        return build_tf_model(lognormal_factory, [-2.0, 0.5])
    if name == "tf_truncated_lognormal":
        return build_tf_model(truncated_lognormal_factory, [-1.0, 0.5, 0.1])
    if name == "tf_gamma":
        return build_tf_model(gamma_factory, [1.3, 3.0])

    raise Exception(f"Model {name} not known.")
