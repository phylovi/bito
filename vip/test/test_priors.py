import numpy as np
from pytest import approx
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
import vip.priors as priors

tf.enable_v2_behavior()
tfd = tfp.distributions

SAMPLE = np.array([[1.0, 2.0, 3.0], [0.26097, 0.0286401, 0.113843]])


def test_log_exp_prior():
    theirs = np.sum(tfp.distributions.Exponential(10).log_prob(SAMPLE).numpy(), axis=1)
    ours = priors.log_exp_prior(SAMPLE)
    assert ours == approx(theirs)
