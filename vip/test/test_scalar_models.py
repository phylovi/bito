import numpy as np
from pytest import approx
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
import vip.scalar_model as models

tf.enable_v2_behavior()
tfd = tfp.distributions


def test_lognormal_log_prob():
    sample = np.array([0.3, 3.14159, 0.0001])
    params = np.array([[0.4, 1.3], [-2.0, 4.0], [3.0, 0.2]])
    which_variables = np.arange(len(sample))
    log_normal = models.LogNormalModel(np.array([0.0, 1.0]), len(sample))
    log_normal.q_params[:, :] = params
    tf_log_normal = models.of_name("tf_lognormal", variable_count=len(sample))
    tf_log_normal.q_params[:, :] = params
    ours = log_normal.log_prob(sample, which_variables=which_variables)
    theirs = tf_log_normal.log_prob(sample, which_variables=which_variables)
    assert ours == approx(theirs)
