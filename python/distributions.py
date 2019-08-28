import numpy as np


class Normal(object):
    def __init__(self, dim):
        self.d = dim

    def sample(self, mu, sigma, n_particles=1):
        if n_particles > 1:
            return np.random.normal(loc=mu,
                                    scale=np.abs(sigma),
                                    size=(n_particles, self.d))
        else:
            return np.random.normal(loc=mu, scale=np.abs(sigma))

    def log_prob(self, x, mu, sigma, grad=False):
        ratio = np.exp(np.log((x - mu)**2) - np.log(sigma**2))
        if grad:
            return ratio / (mu - x)
        else:
            return np.sum(
                -0.5 * np.log(2 * np.pi) - 0.5 * np.log(sigma**2) -
                0.5 * ratio,
                axis=1,
            )

    def reparam_grad(self, x, mu, sigma):
        """The gradient component from the reparametrization trick. In the
        reparametrization trick for the normal distribution, our normal
        variate is mu + sigma * epsilon. Here we take the derivative of that
        with respect to mu and sigma.
        """
        # std_x is epsilon. We translate from x to the epsilon that would have
        # made it. std_x is epsilon. We translate from x to the epsilon that
        # would have made it.
        std_x = (x - mu) / sigma
        return np.ones(mu.shape), std_x

    def lp_param_grad(self, x, mu, sigma):
        """The derivative of a the log of distribution with respect to its
        parameters."""
        ratio = np.exp(np.log((x - mu)**2) - np.log(sigma**2))
        return ratio / (x - mu), -1.0 / sigma * (1.0 - ratio)
