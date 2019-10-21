import numpy as np
from pytest import approx
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp
import vip.scalar_model as models

tf.enable_v2_behavior()
tfd = tfp.distributions

params = np.array([[0.4, 1.3], [-2.0, 4.0], [3.0, 0.2]])
variable_count = params.shape[0]


def test_lognormal_log_prob():
    sample = np.array([0.3, 3.14159, 0.0001])
    which_variables = np.arange(variable_count)
    log_normal = models.LogNormalModel(np.array([0.0, 1.0]), len(sample))
    log_normal.q_params[:, :] = params
    tf_log_normal = models.of_name("tf_lognormal", variable_count=len(sample))
    tf_log_normal.q_params[:, :] = params
    ours = log_normal.log_prob(sample, which_variables=which_variables)
    theirs = tf_log_normal.log_prob(sample, which_variables=which_variables)
    assert ours == approx(theirs)


def test_lognormal_gradients():
    tf.random.set_seed(1)
    particle_count = 8
    tf_log_normal = models.of_name("tf_lognormal", variable_count=variable_count)
    tf_log_normal.q_params[:, :] = params
    px_which_variables = np.array(
        [np.arange(variable_count) for _ in range(particle_count)]
    )
    theirs = tf_log_normal.sample_and_gradients(particle_count, px_which_variables)
    sample = theirs[0]
    log_normal = models.LogNormalModel(np.array([0.0, 1.0]), variable_count)
    log_normal.q_params[:, :] = params
    ours = log_normal.sample_and_gradients(
        particle_count, px_which_variables, prebaked_sample=sample
    )
    for (our_item, their_item) in zip(ours, theirs):
        assert our_item == approx(their_item, rel=1e-5)
