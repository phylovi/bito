import click
import numpy as np

import libsbn
import vip.optimizers
import vip.priors
import vip.scalar_models


class Burrito:
    """A class to wrap an instance and relevant model data.

    Some terminology:

    * We sample "particles" in order to calculate the ELBO and to do stochastic gradient
    ascent. This is a handy terminology because it can be used for trees, or scalar
    parameters, etc. We use the prefix `px_` to designate that the first dimension of
    the given object is across particles. For example, `px_branch_lengths` is a list
    of branch length vectors, where the list is across particles.
    """

    def __init__(
        self,
        *,
        mcmc_nexus_path,
        fasta_path,
        model_name,
        optimizer_name,
        particle_count,
        thread_count=1
    ):
        self.particle_count = particle_count
        self.inst = libsbn.instance("burrito")

        # Read MCMC run to get tree structure.
        self.inst.read_nexus_file(mcmc_nexus_path)
        self.inst.process_loaded_trees()

        # Set up tree likelihood calculation.
        self.inst.read_fasta_file(fasta_path)
        self.inst.make_beagle_instances(thread_count)
        sbn_model = vip.sbn_model.SBNModel(self.inst)
        (branch_lengths, branch_to_split) = self.sample_topology()
        scalar_model = vip.scalar_models.of_name(
            model_name, variable_count=len(branch_lengths)
        )
        self.opt = vip.optimizers.of_name(optimizer_name, sbn_model, scalar_model)

        # PSP: Make choices about branch representations here, eventually.
        self.branch_representations = self.split_based_representations
        # PSP: We will also want to make a choice about using an alternative to
        # split_based_scalar_grad.

    # We want the models to be part of the optimizer because they have state that's
    # closely connected with optimizer state. But we add these convenience functions:

    @property
    def sbn_model(self):
        return self.opt.sbn_model

    @property
    def scalar_model(self):
        return self.opt.scalar_model

    def split_based_representations(self):
        """The ith entry of this array gives the index of the split
        corresponding to the ith branch."""
        return [
            np.array(representation[0])
            for representation in self.inst.get_psp_indexer_representations()
        ]

    def sample_topology(self):
        """Sample a tree, then set up branch length vector and the translation
        from splits to branches and back again."""
        self.inst.sample_trees(1)
        tree = self.inst.tree_collection.trees[0]
        branch_lengths_extended = np.array(tree.branch_lengths, copy=False)
        # Here we are getting a slice that excludes the last (fake) element.
        # Thus we can just deal with the actual branch lengths.
        branch_lengths = branch_lengths_extended[:-1]
        # Now we need to set things up to translate between split indexing and branch
        # indexing.
        # The ith entry of this array gives the index of the split
        # corresponding to the ith branch.
        branch_to_split = np.array(self.inst.get_psp_indexer_representations()[0][0])
        return (branch_lengths, branch_to_split)

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

    def split_based_scalar_grad(
        self, theta_sample, branch_gradients, px_branch_to_split, dg_dpsi, dlog_qg_dpsi
    ):
        """Do a gradient for the scalar parameters in terms of splits.

        See the tex for details. Comments of the form eq:XXX refer to
        equations in the tex.
        See class-level docstring for information about `px_`.
        """
        # Calculate the gradient of the log unnormalized posterior.
        dlogp_dtheta = np.zeros_like(theta_sample)
        for particle_idx, (_, log_grad_raw) in enumerate(branch_gradients):
            # This :-2 is because of the two trailing zeroes that appear at the end of
            # the gradient.
            dlogp_dtheta[particle_idx, :] = np.array(log_grad_raw, copy=False)[:-2]
        dlogp_dtheta += vip.priors.grad_log_exp_prior(theta_sample)
        # Now build up the complete gradient with respect to the variables.
        grad = np.zeros(
            (self.scalar_model.variable_count, self.scalar_model.param_count)
        )
        for particle_idx, branch_to_split in enumerate(px_branch_to_split):
            for branch_idx, variable_idx in enumerate(branch_to_split):
                grad[variable_idx, :] += (
                    # eq:dLdPsi
                    dlogp_dtheta[particle_idx, branch_idx]
                    * dg_dpsi[particle_idx, variable_idx, :]
                    - dlog_qg_dpsi[variable_idx, :]
                )
        return grad

    def gradient_step(self):
        """Take a gradient step."""
        px_branch_lengths = self.sample_topologies(self.particle_count)
        px_branch_representation = self.branch_representations()
        # Sample continuous variables based on the branch representations.
        # PSPs: we'll need to take these branch representations and turn them into
        # the list of variables that they employ.
        (theta_sample, dg_dpsi, dlog_qg_dpsi) = self.scalar_model.sample_and_gradients(
            self.particle_count, px_branch_representation
        )
        # Put the branch lengths in the libsbn instance trees, and get branch gradients.
        for particle_idx, branch_lengths in enumerate(px_branch_lengths):
            branch_lengths[:] = theta_sample[particle_idx, :]
        branch_gradients = self.inst.branch_gradients()
        scalar_grad = self.split_based_scalar_grad(
            theta_sample,
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

    def elbo_estimate(self, particle_count):
        """A naive Monte Carlo estimate of the ELBO."""
        px_branch_lengths = self.sample_topologies(particle_count)
        px_branch_to_split = self.branch_representations()
        # Sample continuous variables based on the branch representations.
        theta_sample = self.scalar_model.sample(particle_count, px_branch_to_split)
        # Put the branch lengths in the libsbn instance trees, and get branch gradients.
        for particle_idx, branch_lengths in enumerate(px_branch_lengths):
            branch_lengths[:] = theta_sample[particle_idx, :]
        px_phylo_log_like = np.array(self.inst.log_likelihoods(), copy=False)
        px_log_prior = vip.priors.log_exp_prior(theta_sample)
        elbo_total = np.sum(px_phylo_log_like + px_log_prior)
        for particle_idx, branch_to_split in enumerate(px_branch_to_split):
            elbo_total -= self.scalar_model.log_prob(
                theta_sample[particle_idx, :], which_variables=branch_to_split
            )
        return elbo_total / particle_count
