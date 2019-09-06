import numpy as np
import scipy as sp


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


def like_weights(q_distribution, phylo_log_likelihoods, x, loc, shape, clip):
    """The weights \tilde w used in the gradient calculation.

    See just before (6) in the 2018 ICLR paper.
    """
    log_prob_ratio = phylo_log_likelihoods - q_distribution.log_prob(x, loc, shape)
    # if clip:
    #     log_prob_ratio = np.clip(log_prob_ratio, -clip, clip)
    return sp.special.softmax(log_prob_ratio)


def param_grad(q_distribution, weights, phylo_gradient, x, loc, shape):
    """The gradient with respect to the parameters of the variational
    distribution using the reparametrization trick.

    We expect x to be two dimensional, with samples on axis 0 and branches on axis 1.

    See (7) of the 2018 ICLR paper.
    """
    assert phylo_gradient.ndim == 2
    # Recall that the branch lengths are parameterized in terms of the variational
    # parameters (epsilon is considered fixed) so when we take their derivative we need
    # to take the branch length derivative:
    d_log_prob_ratio = phylo_gradient - q_distribution.log_prob_grad(x, loc, shape)
    # ... and then the reparametrization derivative due to the chain rule.
    d_reparam_loc, d_reparam_shape = q_distribution.reparam_grad(x, loc, shape)
    # This is the derivative with respect to the variational parametrization itself.
    d_q_loc, d_q_shape = q_distribution.log_prob_param_grad(x, loc, shape)
    # Here are the full derivatives, where the first term is the derivative WRT
    # parameters in the reparameterization trick and the second is WRT the actual
    # variational parameterization.
    d_loc = d_log_prob_ratio * d_reparam_loc - d_q_loc
    d_shape = d_log_prob_ratio * d_reparam_shape - d_q_shape
    # Our multiple samples are laid out on axis 0, so this multiplication on the left
    # reweights them.
    #weights = np.ones(weights.shape) / len(weights)
    return np.matmul(weights, d_loc), np.matmul(weights, d_shape)
