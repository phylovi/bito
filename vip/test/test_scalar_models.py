import numpy as np
from pytest import approx
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
import vip.scalar_models as models

tf.enable_v2_behavior()
tfd = tfp.distributions


def test_lognormal_log_prob():
    x = 0.3
    loc = 0.4
    scale = 1.3
    log_normal = models.LogNormalModel(np.array([loc, scale]), 1)
    ours = log_normal.log_prob(np.array([x]), which_variables=np.array([0]))
    params = tf.constant([x, loc, scale])
    distribution = tfp.distributions.LogNormal(loc=params[1], scale=params[2])
    y = distribution.log_prob(params[0])
    theirs = y.numpy()
    assert ours == approx(theirs)
