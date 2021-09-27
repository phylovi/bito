import numpy as np
from pytest import approx
import bito
import vip.burrito


def test_elbo_innards():
    """ From Mathieu:
    tree.mars.distance mu: -1.728809 sigma: 0.459529 sample: 0.184472
    tree.saturn.distance mu: -2.410943 sigma: 0.748569 sample: 0.027993
    tree.jupiter.distance mu: -2.410977 sigma: 0.748571 sample: 0.045583
    like: -81.446550 prior: 4.327275 logQ: 5.330697
    elbo: -82.449972
    """

    phylo_model_specification = bito.PhyloModelSpecification(
        substitution="JC69", site="constant", clock="strict"
    )

    burro = vip.burrito.Burrito(
        mcmc_nexus_path="data/hello_out.t",
        burn_in_fraction=0,
        fasta_path="data/hello.fasta",
        phylo_model_specification=phylo_model_specification,
        branch_model_name="split",
        scalar_model_name="lognormal",
        optimizer_name="simple",
        particle_count=1,
        thread_count=1,
    )
    branch_model = burro.branch_model

    px_branch_lengths = burro.sample_topologies(1)
    branch_lengths = np.array(px_branch_lengths[0], copy=False)
    theta_sample = np.array([0.184472, 0.027993, 0.045583])
    branch_lengths[:] = theta_sample
    px_theta_sample = np.array([theta_sample])

    mathieu_q_params = np.array(
        [[-1.728809, 0.459529], [-2.410943, 0.748569], [-2.410977, 0.748571]]
    )
    px_branch_representation = branch_model.px_branch_representation()
    branch_rep = px_branch_representation[0]
    # So if the 0th entry of branch_rep is 1, then we are setting the 1th entry of our
    # parameters to the 0th entry of mathieu's, which is in terms of branches.
    branch_model.scalar_model.q_params[branch_rep, :] = mathieu_q_params

    assert np.array(burro.inst.log_likelihoods())[0] == approx(-81.446550)
    assert burro.branch_model.log_prior(px_theta_sample)[0] == approx(4.327275)
    assert burro.branch_model.log_prob(
        px_theta_sample, px_branch_representation
    ) == approx(5.330697, rel=1e-5)
