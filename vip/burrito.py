"""The Burrito class wraps an instance and relevant model data."""
import click
import numpy as np
from scipy.special import logsumexp

import bito
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
        use_vimco=True
    ):
        self.particle_count = particle_count
        self.use_vimco = use_vimco
        self.inst = bito.unrooted_instance("burrito")

        # Read MCMC run to get tree structure.
        self.inst.read_nexus_file(mcmc_nexus_path)
        burn_in_count = int(burn_in_fraction * self.inst.tree_count())
        self.inst.tree_collection.erase(0, burn_in_count)
        self.inst.process_loaded_trees()

        # Set up tree likelihood calculation.
        self.inst.read_fasta_file(fasta_path)
        self.inst.prepare_for_phylo_likelihood(
            phylo_model_specification, thread_count, [], True, particle_count
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
        self.elbo_trace = []

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

    def gradient_step(self, beta_t=1.0):
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
        # Put the branch lengths in the bito instance trees, and get branch gradients.
        for particle_idx, branch_lengths in enumerate(px_branch_lengths):
            branch_lengths[:] = px_theta_sample[particle_idx, :]
        phylo_gradients = self.inst.phylo_gradients()
        scalar_grad = self.branch_model.scalar_grad(
            px_theta_sample,
            phylo_gradients,
            px_branch_representation,
            dg_dpsi,
            dlog_qg_dpsi,
        )
        px_phylo_log_like = np.array(
            [gradient.log_likelihood for gradient in phylo_gradients]
        )
        px_phylo_log_like = beta_t * px_phylo_log_like
        px_log_f = self.px_log_f(
            px_phylo_log_like, px_theta_sample, px_branch_representation
        )
        # Get topology gradients.
        sbn_grad = self.inst.topology_gradients(px_log_f, self.use_vimco)
        self.opt.gradient_step({"scalar_params": scalar_grad, "sbn_params": sbn_grad})

    def gradient_steps(self, step_count):
        betas = np.arange(1, step_count + 1, dtype=np.float)
        betas = np.maximum(betas / step_count, 0.001)
        with click.progressbar(range(step_count), label="Gradient descent") as bar:
            for step in bar:
                self.gradient_step(betas[step])
                self.elbo_trace.append(self.estimate_elbo(self.particle_count))

    def estimate_elbo(self, particle_count):
        """Sample particle_count particles and then make a naive Monte Carlo
        estimate of the ELBO."""
        px_branch_lengths = self.sample_topologies(particle_count)
        px_branch_representation = self.branch_model.px_branch_representation()
        # Sample continuous variables based on the branch representations.
        px_theta_sample = self.branch_model.sample(px_branch_representation)
        # Put the branch lengths in the bito instance trees, and get branch gradients.
        for particle_idx, branch_lengths in enumerate(px_branch_lengths):
            branch_lengths[:] = px_theta_sample[particle_idx, :]
        self.inst.resize_phylo_model_params()
        px_phylo_log_like = np.array(self.inst.log_likelihoods(), copy=False)
        return self.elbo_of_sample(
            px_phylo_log_like,
            px_theta_sample,
            px_branch_representation,
        )

    def elbo_of_sample(
        self, px_phylo_log_like, px_theta_sample, px_branch_representation
    ):
        """A naive Monte Carlo estimate of the ELBO given a sample."""
        px_log_prior = self.branch_model.log_prior(px_theta_sample)
        elbo_total = (
            np.sum(px_phylo_log_like + px_log_prior)
            - np.sum(np.log(self.inst.calculate_sbn_probabilities()))
            - self.branch_model.log_prob(px_theta_sample, px_branch_representation)
        )
        return elbo_total / self.inst.tree_count()

    def px_log_f(self, px_phylo_log_like, px_theta_sample, px_branch_representation):
        """The vector of f values for each tree."""
        px_log_prior = self.branch_model.log_prior(px_theta_sample)
        px_log_sbn_prob = np.log(self.inst.calculate_sbn_probabilities())
        px_branch_log_prob = np.array(
            list(
                self.branch_model.log_prob_generator(
                    px_theta_sample, px_branch_representation
                )
            )
        )
        return px_phylo_log_like + px_log_prior - px_log_sbn_prob - px_branch_log_prob

    def marginal_likelihood_estimate(self, particle_count):
        px_branch_lengths = self.sample_topologies(particle_count)
        px_branch_representation = self.branch_model.px_branch_representation()
        # Sample continuous variables based on the branch representations.
        px_theta_sample = self.branch_model.sample(px_branch_representation)
        # Put the branch lengths in the bito instance trees, and get branch gradients.
        for particle_idx, branch_lengths in enumerate(px_branch_lengths):
            branch_lengths[:] = px_theta_sample[particle_idx, :]
        self.inst.resize_phylo_model_params()
        px_phylo_log_like = np.array(self.inst.log_likelihoods(), copy=False)

        px_log_f = self.px_log_f(
            px_phylo_log_like, px_theta_sample, px_branch_representation
        )

        return logsumexp(px_log_f) - np.log(particle_count)
