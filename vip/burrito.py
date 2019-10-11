import click
import numpy as np

import libsbn
import vip.scalar_models
import vip.optimizers


class Burrito:
    """A class to wrap an instance and relevant model data.

    The current division of labor is that the optimizer handles
    everything after we have sampled a topology, while the burrito can
    sample topologies and then ask the optimizer to update model
    parameters accordingly.
    """

    def __init__(
        self,
        *,
        mcmc_nexus_path,
        fasta_path,
        model_name,
        optimizer_name,
        step_count,
        particle_count
    ):
        self.inst = libsbn.instance("burrito")

        # Read MCMC run to get tree structure.
        self.inst.read_nexus_file(mcmc_nexus_path)
        self.inst.process_loaded_trees()

        # Set up tree likelihood calculation.
        self.inst.read_fasta_file(fasta_path)
        self.inst.make_beagle_instances(1)
        # It's important to do tree sampling here so that self.branch_lengths
        # etc gets set up.
        self.sample_topology()

        sbn_model = vip.sbn_model.SBNModel(self.inst)
        scalar_model = vip.scalar_models.of_name(
            model_name,
            # ONEBUG
            variable_count=1 + len(self.branch_lengths),
            particle_count=particle_count,
        )
        self.opt = vip.optimizers.of_name(optimizer_name, sbn_model, scalar_model)

    @staticmethod
    def log_exp_prior(x, rate=10):
        return np.log(rate) * x.shape[1] - rate * np.sum(x, axis=1)

    @staticmethod
    def grad_log_exp_prior(x, rate=10):
        return -rate

    def sample_topology(self):
        """Sample a tree, then set up branch length vector and the translation
        from splits to branches and back again."""
        self.inst.sample_trees(1)
        tree = self.inst.tree_collection.trees[0]
        branch_lengths_extended = np.array(tree.branch_lengths, copy=False)
        # Here we are getting a slice that excludes the last (fake) element.
        # Thus we can just deal with the actual branch lengths.
        self.branch_lengths = branch_lengths_extended[:-1]

        # Note: in the following arrays, the particles are laid out along axis 0 and the
        # splits are laid out along axis 1.
        # Now we need to set things up to translate between split indexing and branch
        # indexing.
        # The ith entry of this array gives the index of the split
        # corresponding to the ith branch.
        self.branch_to_split = np.array(
            self.inst.get_psp_indexer_representations()[0][0]
        )

    def log_like_with(self, branch_lengths):
        self.branch_lengths[:] = branch_lengths
        return np.array(self.inst.log_likelihoods())[0]

    def grad_log_like_with(self, branch_lengths):
        self.branch_lengths[:] = branch_lengths
        _, log_grad = self.inst.branch_gradients()[0]
        # This :-2 is because of the two trailing zeroes that appear at the end of
        # the gradient.
        return np.array(log_grad)[:-2]

    def phylo_log_like(self, branch_lengths_arr):
        """Calculate phylogenetic log likelihood for each of the split length
        assignments laid out along axis 1."""
        return np.apply_along_axis(self.log_like_with, 1, branch_lengths_arr)

    def grad_phylo_log_like(self, branch_lengths_arr):
        return np.apply_along_axis(self.grad_log_like_with, 1, branch_lengths_arr)

    def phylo_log_upost(self, branch_lengths_arr):
        """The unnormalized phylogenetic posterior with an Exp(10) prior."""
        return self.phylo_log_like(branch_lengths_arr) + Burrito.log_exp_prior(
            branch_lengths_arr
        )

    def grad_phylo_log_upost(self, branch_lengths_arr):
        """The unnormalized phylogenetic posterior with an Exp(10) prior."""
        return self.grad_phylo_log_like(
            branch_lengths_arr
        ) + Burrito.grad_log_exp_prior(branch_lengths_arr)

    def gradient_steps(self, step_count):
        with click.progressbar(range(step_count), label="Gradient descent") as bar:
            for step in bar:
                # SOON: self.sample_topology()
                if not self.opt.gradient_step(
                    self.phylo_log_upost,
                    self.grad_phylo_log_upost,
                    self.branch_to_split,
                ):
                    raise Exception("ELBO is not finite. Stopping.")

    def elbo_estimate(self, particle_count=None):
        """A naive Monte Carlo estimate of the ELBO."""
        if particle_count is None:
            particle_count = self.particle_count
        z, log_prob = self.opt.scalar_model.sample(
            particle_count, which_variables=self.branch_to_split, log_prob=True
        )
        return (np.sum(self.phylo_log_upost(z) - log_prob)) / particle_count
