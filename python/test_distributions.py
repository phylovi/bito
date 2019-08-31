import distributions

from pytest import approx
import tensorflow as tf
import tensorflow_probability as tfp

tf.enable_eager_execution()

# Test for Gamma
x = 0.3
alpha = 0.55
beta = 0.6
params = tf.constant([x, alpha, beta])
with tf.GradientTape() as g:
    g.watch(params)
    gamma = tfp.distributions.Gamma(concentration=params[1], rate=params[2])
    y = gamma.log_prob(params[0])
gradient = g.gradient(y, params).numpy()
g = distributions.Gamma(1)
assert gradient[0] == approx(g.log_prob_grad(x, alpha, beta))

# Test for Normal
x = 0.3
loc = 0.55
shape = 0.6
params = tf.constant([x, loc, shape])
with tf.GradientTape() as g:
    g.watch(params)
    normal = tfp.distributions.Normal(loc=params[1], scale=params[2])
    y = normal.log_prob(params[0])
gradient = g.gradient(y, params).numpy()
d = distributions.Normal(1)
assert gradient[0] == approx(d.log_prob_grad(x, loc, shape))
assert (gradient[1], gradient[2]) == approx(d.log_prob_param_grad(x, loc, shape))
