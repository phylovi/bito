import click
import numpy as np

import libsbn
import vip.branch_model
import vip.optimizers
import vip.priors


class Burrito:
    """A class to wrap an instance and relevant model data.

    Some terminology:

    * We sample "particles" in order to calculate the ELBO and to do stochastic gradient
    ascent. This is a handy terminology because it can be used for trees, or scalar
    parameters, etc. We use the prefix `px_` to designate that the first dimension of
    the given object is across particles. For example, `px_branch_lengths` is a list
    of branch length vectors, where the list is across particles.

    * particle_count is the particle count to be used for gradient calculation
    """

    def __init__(
        self,
        *,
        mcmc_nexus_path,
        burn_in_fraction,
        fasta_path,
        phylo_model_specification,
        branch_model_name,
        scalar_model_name,
        optimizer_name,
        particle_count,
        thread_count=1,
    ):
        self.particle_count = particle_count
        self.inst = libsbn.instance("burrito")

        # Read MCMC run to get tree structure.
        self.inst.read_nexus_file(mcmc_nexus_path)
        burn_in_count = int(burn_in_fraction * self.inst.tree_count())
        self.inst.tree_collection.erase(0, burn_in_count)
        self.inst.process_loaded_trees()

        # Set up tree likelihood calculation.
        self.inst.read_fasta_file(fasta_path)
        self.inst.prepare_for_phylo_likelihood(
            phylo_model_specification, thread_count, particle_count
        )
        sbn_model = vip.sbn_model.SBNModel(self.inst)
        self.branch_model = vip.branch_model.of_name(
            branch_model_name, scalar_model_name, self.inst
        )
        self.opt = vip.optimizers.of_name(
            optimizer_name,
            sbn_model,
            self.branch_model.scalar_model,
            self.estimate_elbo,
        )

    # We want the models to be part of the optimizer because they have state that's
    # closely connected with optimizer state. But we add these convenience functions:

    @property
    def sbn_model(self):
        return self.opt.sbn_model

    def sample_topologies(self, count):
        """Sample trees into the instance and return the np'd version of their
        branch length vectors."""
        self.inst.sample_trees(count)
        # Here we are getting a slice that excludes the last (fake) element.
        # Thus we can just deal with the actual branch lengths.
        return [
            np.array(tree.branch_lengths, copy=False)[:-1]
            for tree in self.inst.tree_collection.trees
        ]

    def gradient_step(self):
        """Take a gradient step."""
        px_branch_lengths = self.sample_topologies(self.particle_count)
        px_branch_representation = self.branch_model.px_branch_representation()
        # This design may seem a little strange, in that we separate out the
        # branch_model part of the gradients and then consume them via branch_model
        # later. This is driven by wanting to support scalar models (such as the TF
        # models) that require sampling and gradient calculation at the same time.
        (
            px_theta_sample,
            dg_dpsi,
            dlog_qg_dpsi,
        ) = self.branch_model.sample_and_gradients(px_branch_representation)
        # Put the branch lengths in the libsbn instance trees, and get branch gradients.
        for particle_idx, branch_lengths in enumerate(px_branch_lengths):
            branch_lengths[:] = px_theta_sample[particle_idx, :]
        branch_gradients = self.inst.branch_gradients()
        scalar_grad = self.branch_model.scalar_grad(
            px_theta_sample,
            branch_gradients,
            px_branch_representation,
            dg_dpsi,
            dlog_qg_dpsi,
        )
        self.opt.gradient_step(scalar_grad)

    def gradient_steps(self, step_count):
        with click.progressbar(range(step_count), label="Gradient descent") as bar:
            for step in bar:
                self.gradient_step()

    def estimate_elbo(self, particle_count):
        """A naive Monte Carlo estimate of the ELBO."""
        px_branch_lengths = self.sample_topologies(particle_count)
        px_branch_representation = self.branch_model.px_branch_representation()
        # Sample continuous variables based on the branch representations.
        px_theta_sample = self.branch_model.sample(px_branch_representation)
        # Put the branch lengths in the libsbn instance trees, and get branch gradients.
        for particle_idx, branch_lengths in enumerate(px_branch_lengths):
            branch_lengths[:] = px_theta_sample[particle_idx, :]
        self.inst.resize_phylo_model_params()
        px_phylo_log_like = np.array(self.inst.log_likelihoods(), copy=False)
        px_log_prior = self.branch_model.log_prior(px_theta_sample)
        elbo_total = np.sum(
            px_phylo_log_like + px_log_prior
        ) - self.branch_model.log_prob(px_theta_sample, px_branch_representation)
        return elbo_total / particle_count
