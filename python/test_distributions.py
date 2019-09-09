import branch_length
import distributions

import numpy as np
from pytest import approx
import tensorflow as tf
import tensorflow_probability as tfp

tf.enable_eager_execution()

def detensor_pair(pair):
    return pair[0][0,0], pair[1][0,0]

def test_distribution_gradients():
    x = 0.3
    x_arr = np.array([[x]])

    # Test for Gamma
    alpha = np.array([0.55])
    beta = np.array([0.6])
    params = tf.constant([x, alpha[0], beta[0]])
    with tf.GradientTape() as g:
        g.watch(params)
        tf_distribution = tfp.distributions.Gamma(
            concentration=params[1], rate=params[2]
        )
        y = tf_distribution.log_prob(params[0])
    gradient = g.gradient(y, params).numpy()
    g = distributions.Gamma(1)
    assert y.numpy() == approx(g.log_prob(x_arr, alpha, beta)[0])
    assert gradient[0] == approx(g.log_prob_grad(x_arr, alpha, beta))

    # Test for Normal
    loc = np.array([0.55])
    shape = np.array([0.6])
    params = tf.constant([x, loc[0], shape[0]])
    with tf.GradientTape() as g:
        g.watch(params)
        tf_distribution = tfp.distributions.Normal(loc=params[1], scale=params[2])
        y = tf_distribution.log_prob(params[0])
    gradient = g.gradient(y, params).numpy()
    d = distributions.Normal(1)
    assert y.numpy() == approx(d.log_prob(x_arr, loc, shape)[0])
    assert gradient[0] == approx(d.log_prob_grad(x_arr, loc, shape))
    assert (gradient[1], gradient[2]) == approx(detensor_pair(d.log_prob_param_grad(x_arr, loc, shape)))

    # Test for LogNormal
    loc = np.array([0.55])
    shape = np.array([0.6])
    params = tf.constant([x, loc[0], shape[0]])
    with tf.GradientTape() as g:
        g.watch(params)
        tf_distribution = tfp.distributions.LogNormal(loc=params[1], scale=params[2])
        y = tf_distribution.log_prob(params[0])
    gradient = g.gradient(y, params).numpy()
    d = distributions.LogNormal(1)
    assert y.numpy() == approx(d.log_prob(x_arr, loc, shape)[0])
    assert gradient[0] == approx(d.log_prob_grad(x_arr, loc, shape))
    assert (gradient[1], gradient[2]) == approx(detensor_pair(d.log_prob_param_grad(x_arr, loc, shape)))


def test_log_ratio_gradient():
    """
    Here we test our like_weights and param_grad using TF.
    """
    true_loc_val = -2.0
    true_shape_val = 1.4
    true_normal = tfp.distributions.Normal(loc=true_loc_val, scale=true_shape_val)
    # A variety of epsilons.
    epsilon = tf.constant([-0.1, 0.14, -0.51, -2.0, 1.8])
    with tf.GradientTape() as g:
        mu = tf.constant(1.1)
        sigma = tf.constant(1.0)
        g.watch(mu)
        g.watch(sigma)
        tf_x = mu + sigma * epsilon
        # This is the log of the full sum of ratios as in the equation just before (7)
        # in the 2018 ICLR paper.
        y = tf.math.log(
            # In principle we should have a 1/K term here, but it disappears in the log
            # grad.
            tf.math.reduce_sum(
                true_normal.prob(tf_x)
                / tfp.distributions.Normal(loc=mu, scale=sigma).prob(tf_x)
            )
        )
        x_arr = np.array([tf_x.numpy()]).transpose()
        tf_gradient = [grad.numpy() for grad in g.gradient(y, [mu, sigma])]
    # Distribution we wish to approximate.
    d = distributions.Normal(1)
    true_loc = np.array([true_loc_val])
    true_shape = np.array([true_shape_val])
    # Variational distribution
    q = distributions.Normal(1)
    loc = np.array([1.1])
    shape = np.array([1.0])
    # Here's the gradient as per (7) in the 2018 ICLR paper.
    def complete_grad(x, loc, shape):
        weights = branch_length.like_weights(
            q, d.log_prob(x, true_loc, true_shape), x, loc, shape, None
        )
        return branch_length.param_grad(
            q, weights, d.log_prob_grad(x, true_loc, true_shape), x, loc, shape
        )
    loc_grad, shape_grad = complete_grad(x_arr, loc, shape)
    assert tf_gradient[0] == approx(loc_grad, rel=1e-5)
    assert tf_gradient[1] == approx(shape_grad, rel=1e-5)
