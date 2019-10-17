import click
import numpy as np

import libsbn
import vip.optimizers
import vip.priors
import vip.scalar_models


class Burrito:
    """A class to wrap an instance and relevant model data.

    TODO overview.
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
        self.sbn_model = vip.sbn_model.SBNModel(self.inst)
        (branch_lengths, branch_to_split) = self.sample_topology()
        scalar_model = vip.scalar_models.of_name(
            model_name,
            variable_count=len(branch_lengths),
            particle_count=particle_count,
        )
        self.opt = vip.optimizers.of_name(optimizer_name, self.sbn_model, scalar_model)
        self.scalar_model = self.opt.scalar_model

        # Make choices about branch representations here, eventually.
        self.branch_representations = self.split_based_representations

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

    def sample_topologies(self):
        """Sample trees into the instance and return the np'd version of their
        branch length vectors."""
        self.inst.sample_trees(self.particle_count)
        # Here we are getting a slice that excludes the last (fake) element.
        # Thus we can just deal with the actual branch lengths.
        # TODO explain what arr means, and/or change it.
        return [
            np.array(tree.branch_lengths, copy=False)[:-1]
            for tree in self.inst.tree_collection.trees
        ]

    def gradient_step(self):
        """Take a gradient step."""
        branch_lengths_arr = self.sample_topologies()
        branch_to_split_arr = self.branch_representations()
        # For PSPs, we'll need to take these branch representations and turn them into
        # the list of variables that they employ.
        (
            theta_sample,
            dg_dpsi,
            dlog_sum_q_dpsi,
        ) = self.scalar_model.sample_and_gradients(branch_to_split_arr)
        # Put the branch lengths in the libsbn instance trees.
        for particle_idx, branch_lengths in enumerate(branch_lengths_arr):
            branch_lengths[:] = theta_sample[particle_idx, :]
        gradient_result = self.inst.branch_gradients()
        # Calculate the gradient of the log unnormalized posterior.
        dlogp_dtheta = np.zeros_like(theta_sample)
        for particle_idx, (_, log_grad_raw) in enumerate(gradient_result):
            # This :-2 is because of the two trailing zeroes that appear at the end of
            # the gradient.
            dlogp_dtheta[particle_idx, :] = np.array(log_grad_raw, copy=False)[:-2]
        dlogp_dtheta += vip.priors.grad_log_exp_prior(theta_sample)
        # Now build up the complete gradient with respect to the variables.
        vars_grad = np.zeros(
            (self.scalar_model.variable_count, self.scalar_model.param_count)
        )
        for particle_idx, branch_to_split in enumerate(branch_to_split_arr):
            for branch_index, variable_index in enumerate(branch_to_split):
                for param_index in range(self.scalar_model.param_count):
                    vars_grad[variable_index, param_index] += (
                        dlogp_dtheta[particle_idx, branch_index]
                        * dg_dpsi[particle_idx, variable_index, param_index]
                        # TODO I don't understand why this particle count division.
                        - dlog_sum_q_dpsi[variable_index, param_index]
                        / self.particle_count
                    )
        self.opt.gradient_step(vars_grad)

    def gradient_steps(self, step_count):
        with click.progressbar(range(step_count), label="Gradient descent") as bar:
            for step in bar:
                self.gradient_step()

    def elbo_estimate(self, particle_count=None):
        """A naive Monte Carlo estimate of the ELBO."""
        # TODO update this with sample_topologies
        if particle_count is None:
            particle_count = self.particle_count
        (branch_lengths, branch_to_split) = self.sample_topology()
        theta = self.opt.scalar_model.sample(
            particle_count, which_variables=branch_to_split
        )
        log_prob = self.opt.scalar_model.log_prob(
            theta, which_variables=branch_to_split
        )

        def log_like_with(in_branch_lengths):
            branch_lengths[:] = in_branch_lengths
            return np.array(self.inst.log_likelihoods())[0]

        def phylo_log_upost(branch_lengths_arr):
            """The unnormalized phylogenetic posterior with an Exp(10)
            prior."""
            return np.apply_along_axis(
                log_like_with, 1, branch_lengths_arr
            ) + vip.priors.log_exp_prior(branch_lengths_arr)

        return (np.sum(phylo_log_upost(theta) - log_prob)) / particle_count
