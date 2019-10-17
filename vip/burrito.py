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
        self, *, mcmc_nexus_path, fasta_path, model_name, optimizer_name, particle_count
    ):
        self.particle_count = particle_count
        self.inst = libsbn.instance("burrito")

        # Read MCMC run to get tree structure.
        self.inst.read_nexus_file(mcmc_nexus_path)
        self.inst.process_loaded_trees()

        # Set up tree likelihood calculation.
        self.inst.read_fasta_file(fasta_path)
        self.inst.make_beagle_instances(1)
        self.sbn_model = vip.sbn_model.SBNModel(self.inst)
        # TODO cut this in favor of using an indexer
        (branch_lengths, branch_to_split) = self.sample_topology()
        scalar_model = vip.scalar_models.of_name(
            model_name,
            variable_count=len(branch_lengths),
            particle_count=particle_count,
        )
        self.opt = vip.optimizers.of_name(optimizer_name, self.sbn_model, scalar_model)
        self.scalar_model = self.opt.scalar_model

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
        """Sample a tree, then set up branch length vector and the translation
        from splits to branches and back again."""
        self.inst.sample_trees(self.particle_count)
        # Here we are getting a slice that excludes the last (fake) element.
        # Thus we can just deal with the actual branch lengths.
        branch_lengths_arr = [
            np.array(tree.branch_lengths, copy=False)[:-1]
            for tree in self.inst.tree_collection.trees
        ]
        # The ith entry of this array gives the index of the split
        # corresponding to the ith branch.
        branch_to_split_arr = [
            np.array(representation[0])
            for representation in self.inst.get_psp_indexer_representations()
        ]
        return (branch_lengths_arr, branch_to_split_arr)

    def gradient_step(self):
        """Take a gradient step."""
        (branch_lengths, branch_to_split) = self.sample_topology()
        # Temporarily a constant array.
        which_variables_arr = [branch_to_split for _ in range(self.particle_count)]
        (
            theta_sample,
            dg_dpsi,
            dlog_sum_q_dpsi,
        ) = self.scalar_model.sample_and_gradients(which_variables_arr)
        # Set branch lengths using the scalar model sample

        def grad_log_like_with(in_branch_lengths):
            branch_lengths[:] = in_branch_lengths
            _, log_grad = self.inst.branch_gradients()[0]
            # This :-2 is because of the two trailing zeroes that appear at the end of
            # the gradient.
            return np.array(log_grad)[:-2]

        def grad_phylo_log_upost(branch_lengths_arr):
            return np.apply_along_axis(
                grad_log_like_with, 1, branch_lengths_arr
            ) + vip.priors.grad_log_exp_prior(branch_lengths_arr)

        grad_log_p = grad_phylo_log_upost(theta_sample)
        vars_grad = np.zeros(
            (self.scalar_model.variable_count, self.scalar_model.param_count)
        )
        for branch_index, variable_index in enumerate(branch_to_split):
            for param_index in range(self.scalar_model.param_count):
                vars_grad[variable_index, param_index] += (
                    np.sum(
                        grad_log_p[:, branch_index]
                        * dg_dpsi[:, variable_index, param_index],
                        axis=0,
                    )
                    - dlog_sum_q_dpsi[variable_index, param_index]
                )
        self.opt.gradient_step(vars_grad)

    def gradient_steps(self, step_count):
        with click.progressbar(range(step_count), label="Gradient descent") as bar:
            for step in bar:
                self.gradient_step()

    def elbo_estimate(self, particle_count=None):
        """A naive Monte Carlo estimate of the ELBO."""
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
            """The unnormalized phylogenetic posterior with an Exp(10) prior."""
            return np.apply_along_axis(
                log_like_with, 1, branch_lengths_arr
            ) + vip.priors.log_exp_prior(branch_lengths_arr)

        return (np.sum(phylo_log_upost(theta) - log_prob)) / particle_count
