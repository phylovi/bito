from collections import OrderedDict


def param_grad(q_distribution, phylo_gradient, x, loc, shape):
    """The gradient with respect to the parameters of the variational
    distribution using the reparametrization trick. See (7) of the 2018 ICLR
    paper.

    This explodes if x == loc.
    """
    # Recall that the branch lengths are parameterized in terms of the variational
    # parameters (epsilon is considered fixed) so when we take their derivative we need
    # to take the branch length derivative then do the chain rule. The product of these
    # first two terms are that.
    d_log_probs_ratio = phylo_gradient - q_distribution.log_prob_grad(x, loc, shape)
    d_reparam_loc, d_reparam_shape = q_distribution.reparam_grad(x, loc, shape)
    # This is the derivative with respect to the variational parametrization itself.
    d_q_loc, d_q_shape = q_distribution.log_prob_param_grad(x, loc, shape)
    d_loc = d_log_probs_ratio * d_reparam_loc - d_q_loc
    d_shape = d_log_probs_ratio * d_reparam_shape - d_q_shape

    OrderedDict(
            [
                ("phylo_gradient", phylo_gradient),
                ("q_log_prob_grad", q_distribution.log_prob_grad(x, loc, shape)),
                ("d_log_probs_ratio", d_log_probs_ratio),
                ("d_reparam_loc", d_reparam_loc),
                ("d_reparam_shape", d_reparam_shape),
                ("d_q_loc", d_q_loc),
                ("d_q_shape", d_q_shape),
                ("d_loc", d_loc),
                ("d_shape", d_shape),
            ]
        )

    return d_loc, d_shape
