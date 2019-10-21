import numpy as np
from pytest import approx
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
import vip.scalar_model as models

tf.enable_v2_behavior()
tfd = tfp.distributions


def test_lognormal_log_prob():
    sample = np.array([0.3, 3.14159, 0.0001])
    loc = 0.4
    scale = 1.3
    log_normal = models.LogNormalModel(np.array([loc, scale]), len(sample))
    tf_log_normal = models.of_name("tf_lognormal", variable_count=len(sample))
    tf_log_normal.q_params[:, 0] = loc
    tf_log_normal.q_params[:, 1] = scale
    which_variables = np.arange(len(sample))
    ours = log_normal.log_prob(sample, which_variables=which_variables)
    theirs = tf_log_normal.log_prob(sample, which_variables=which_variables)
    assert ours == approx(theirs)
