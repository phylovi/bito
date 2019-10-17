import numpy as np


def log_exp_prior(x, rate=10):
    # This commented version is the correct one, and will be fixed when I address
    # issue #116.
    # return np.log(rate) * x.shape[1] - rate * np.sum(x, axis=1)
    return np.log(rate) - np.sum(rate * x, axis=1)


def grad_log_exp_prior(x, rate=10):
    return -rate
