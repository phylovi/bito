import numpy as np


def softmax(arrA):
    """Apply softmax function to an array-like.

    Softmax function defined as exp(z_i)/sum_j(exp(z_j)).
    See https://en.wikipedia.org/wiki/Softmax_function

    :param arrA: list or array-like of real numbers.
    :return: array of values in the interval (0, 1) that sum to 1.0.
    """
    max_arrA = np.max(arrA)
    if max_arrA != -np.inf:
        normalizing_constant = np.sum(np.exp(arrA - max_arrA))
        return np.exp(arrA - max_arrA) / normalizing_constant
    else:
        return np.ones(arrA.shape) / len(arrA)


def like_weights(q_distribution, phylo_log_like, x, loc, shape, clip):
    """The gradient with respect to the parameters of the variational
    distribution using the reparametrization trick.

    See (7) of the 2018 ICLR paper.
    """
    # Recall that the branch lengths are parameterized in terms of the variational
    # parameters (epsilon is considered fixed) so when we take their derivative we need
    # to take the branch length derivative then do the chain rule. The product of these
    # first two terms are that.
    log_prob_ratio = phylo_log_like - q_distribution.log_prob(x, loc, shape)
    if clip:
        log_prob_ratio = np.clip(log_prob_ratio, -clip, clip)
    return softmax(log_prob_ratio)


def param_grad(q_distribution, weights, phylo_gradient, x, loc, shape):
    """The gradient with respect to the parameters of the variational
    distribution using the reparametrization trick. See (7) of the 2018 ICLR
    paper.

    TODO: This explodes if x == loc.
    """
    # Recall that the branch lengths are parameterized in terms of the variational
    # parameters (epsilon is considered fixed) so when we take their derivative we need
    # to take the branch length derivative then do the chain rule. The product of these
    # first two terms are that.
    d_log_prob_ratio = phylo_gradient - q_distribution.log_prob_grad(x, loc, shape)
    d_reparam_loc, d_reparam_shape = q_distribution.reparam_grad(x, loc, shape)
    # This is the derivative with respect to the variational parametrization itself.
    d_q_loc, d_q_shape = q_distribution.log_prob_param_grad(x, loc, shape)
    d_loc = d_log_prob_ratio * d_reparam_loc - d_q_loc
    d_shape = d_log_prob_ratio * d_reparam_shape - d_q_shape
    reweight = lambda a: np.dot(np.array([weights]), a)[0]
    return reweight(d_loc), reweight(d_shape)
