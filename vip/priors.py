import numpy as np


def log_exp_prior(theta_sample, rate=10):
    """The log of the Exponential(rate) probability density.

    theta_sample should have particles on axis 0 and variables on axis
    1.
    """
    assert theta_sample.ndim == 2
    return np.log(rate) * theta_sample.shape[1] - rate * np.sum(theta_sample, axis=1)


def grad_log_exp_prior(theta_sample, rate=10):
    """The gradient of the log Exponential(rate) probability density."""
    return -rate
