import numpy as np


class Normal(object):
    """This class implements a multidimensional normal distribution such that
    each dimension has a mu and a sigma.

    The given dimension is laid out on axis 1, which in our application
    is across branches of a tree.

    pdf(x; mu, sigma) = exp(-0.5 (x - mu)**2 / sigma**2) / Z
    Z = (2 pi sigma**2)**0.5
    """

    def __init__(self, dim):
        self.dim = dim

    def sample(self, mu, sigma, n_particles=1):
        """mu and sigma should be dim-dimensional."""
        assert len(mu) == self.dim
        assert len(sigma) == self.dim
        if n_particles > 1:
            return np.random.normal(
                loc=mu, scale=np.abs(sigma), size=(n_particles, self.dim)
            )
        else:
            return np.random.normal(loc=mu, scale=np.abs(sigma))

    def log_prob(self, x, mu, sigma):
        """Given an ensemble of dim-dimensional samples, calculate the joint
        probability (across dim) for each one.

        Thus, x is an array that is samples on axis 0 and branches on
        axis 1.
        """
        assert(x.ndim == 2)
        ratio = self._internal_ratio(x, mu, sigma)
        return np.sum(
            -0.5 * np.log(2 * np.pi) - 0.5 * np.log(sigma ** 2) - 0.5 * ratio, axis=1
        )

    def log_prob_grad(self, x, mu, sigma):
        """The gradient of log_prob."""
        return (x - mu) / sigma ** 2

    def log_prob_param_grad(self, x, mu, sigma):
        """The derivative of a the log of distribution with respect to its
        parameters."""
        ratio = self._internal_ratio(x, mu, sigma)
        return ratio / (x - mu), -1.0 / sigma * (1.0 - ratio)

    def reparam_grad(self, x, mu, sigma):
        """The gradient component from the reparametrization trick. In the
        reparametrization trick for the normal distribution, our normal.

        variate is mu + sigma * epsilon. Here we take the derivative of that
        with respect to mu and sigma.
        """
        # std_x is epsilon. We translate from x to the epsilon that would have
        # made it. std_x is epsilon. We translate from x to the epsilon that
        # would have made it.
        std_x = (x - mu) / sigma
        return np.ones(mu.shape), std_x

    def _internal_ratio(self, x, mu, sigma):
        """(x - mu)**2 / sigma**2

        Cheng used
        ratio = np.exp(np.log((x - mu) ** 2) - np.log(sigma ** 2))
        """
        return (x - mu) ** 2 / sigma ** 2
