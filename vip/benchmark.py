import os

import numpy as np
import pandas as pd
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp

import libsbn
import vip.continuous_parameter_models as models
import vip.optimizers

tf.enable_v2_behavior()
tfd = tfp.distributions


# Documentation is in the CLI.
# `*` forces everything after to be keyword-only.
def fixed(data_path, *, model_name, step_count, particle_count):
    data_path = os.path.normpath(data_path)
    data_id = os.path.basename(data_path)
    mcmc_nexus_path = os.path.join(data_path, data_id + "_out.t")
    fasta_path = os.path.join(data_path, data_id + ".fasta")
    burn_in_fraction = 0.1

    inst = libsbn.instance("charlie")

    # Read MCMC run and get branch lengths.
    inst.read_nexus_file(mcmc_nexus_path)
    inst.process_loaded_trees()
    burn_in_count = int(burn_in_fraction * inst.tree_count())
    mb_split_lengths_np = np.array([np.array(a) for a in inst.split_lengths()])
    mb_split_lengths = pd.DataFrame(mb_split_lengths_np[:, burn_in_count:].transpose())
    last_sampled_split_lengths = mb_split_lengths.iloc[-1].to_numpy()
    mb_split_lengths["total"] = mb_split_lengths.sum(axis=1)

    # Set up tree likelihood calculation.
    inst.read_fasta_file(fasta_path)
    inst.make_beagle_instances(1)
    # This is a "fake" tree sample because there is only one tree in the support of the
    # SBN, but it's the easiest way to get the tree back.
    inst.sample_trees(1)
    tree = inst.tree_collection.trees[0]
    branch_lengths_extended = np.array(tree.branch_lengths, copy=False)
    # Here we are getting a slice that excludes the last (fake) element.
    # Thus we can just deal with the actual branch lengths.
    branch_lengths = branch_lengths_extended[:-1]
    branch_lengths[:] = 0.1

    # Now we need to set things up to translate between split indexing and branch
    # indexing.
    # The ith entry of this array gives the index of the split corresponding to the ith
    # branch.
    branch_to_split = np.array(inst.get_psp_indexer_representations()[0][0])
    # The ith entry of this array gives the index of the branch corresponding to the ith
    # split.
    split_to_branch = np.copy(branch_to_split)
    for branch in range(len(branch_to_split)):
        split_to_branch[branch_to_split[branch]] = branch

    def translate_branches_to_splits(branch_vector):
        """The ith entry of the array returned by this function is the entry of
        branch_vector corresponding to the ith split."""
        return branch_vector[split_to_branch]

    def translate_splits_to_branches(split_vector):
        """The ith entry of the array returned by this function is the entry of
        split_vector corresponding to the ith branch."""
        return split_vector[branch_to_split]

    def log_like_with(split_lengths, grad=False):
        saved_branch_lengths = branch_lengths.copy()
        for branch in range(len(branch_lengths)):
            branch_lengths[branch] = split_lengths[branch_to_split[branch]]
        if grad:
            _, log_grad = inst.branch_gradients()[0]
            result = translate_branches_to_splits(np.array(log_grad)[:-2])
        else:
            result = np.array(inst.log_likelihoods())[0]
            branch_lengths[:] = saved_branch_lengths
        return result

    def phylo_log_like(x_arr):
        """Calculate phylogenetic log likelihood for each of the branch length
        assignments laid out along axis 1."""
        return np.apply_along_axis(log_like_with, 1, x_arr)

    def grad_phylo_log_like(x_arr):
        return np.apply_along_axis(lambda x: log_like_with(x, grad=True), 1, x_arr)

    def log_exp_prior(x, rate=10):
        return np.log(rate) - np.sum(rate * x, axis=1)

    def grad_log_exp_prior(x, rate=10):
        return -rate

    def phylo_log_upost(x_arr):
        """The unnormalized phylogenetic posterior with an Exp(10) prior."""
        return phylo_log_like(x_arr) + log_exp_prior(x_arr)

    def grad_phylo_log_upost(x_arr):
        """The unnormalized phylogenetic posterior with an Exp(10) prior."""
        return grad_phylo_log_like(x_arr) + grad_log_exp_prior(x_arr)

    m = models.of_name(
        model_name, variable_count=len(branch_lengths), particle_count=particle_count
    )
    m.mode_match(last_sampled_split_lengths)
    m.elbo_estimate(phylo_log_upost, particle_count=1000)

    opt = vip.optimizers.AdaptiveStepsizeOptimizer(m)
    opt.gradient_steps(phylo_log_upost, grad_phylo_log_upost, step_count)
    opt_trace = pd.DataFrame({"elbo": opt.trace}).reset_index()

    fit_sample = pd.DataFrame(m.sample(len(mb_split_lengths)))
    fit_sample["total"] = fit_sample.sum(axis=1)
    fit_sample["type"] = "vb"
    mb_split_lengths["type"] = "mcmc"
    fitting_results = pd.concat(
        [fit_sample.melt(id_vars="type"), mb_split_lengths.melt(id_vars="type")]
    )
    fitting_results["variable"] = fitting_results["variable"].astype(str)

    return opt_trace, fitting_results
