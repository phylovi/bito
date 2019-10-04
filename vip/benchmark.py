import os
import timeit

import numpy as np
import pandas as pd
import tensorflow.compat.v2 as tf
import tensorflow_probability as tfp

import libsbn
import vip.continuous_models
import vip.optimizers

tf.enable_v2_behavior()
tfd = tfp.distributions


# Documentation is in the CLI.
# `*` forces everything after to be keyword-only.
def fixed(data_path, *, model_name, optimizer_name, step_count, particle_count):
    data_path = os.path.normpath(data_path)
    data_id = os.path.basename(data_path)
    mcmc_nexus_path = os.path.join(data_path, data_id + "_out.t")
    fasta_path = os.path.join(data_path, data_id + ".fasta")
    burn_in_fraction = 0.1
    particle_count_for_final_elbo_estimate = 2000

    inst = libsbn.instance("charlie")

    # Read MCMC run and get branch lengths.
    inst.read_nexus_file(mcmc_nexus_path)
    inst.process_loaded_trees()
    burn_in_count = int(burn_in_fraction * inst.tree_count())
    mcmc_split_lengths_np = np.array([np.array(a) for a in inst.split_lengths()])
    mcmc_split_lengths = pd.DataFrame(
        mcmc_split_lengths_np[:, burn_in_count:].transpose()
    )
    last_sampled_split_lengths = mcmc_split_lengths.iloc[-1].to_numpy()
    mcmc_split_lengths["total"] = mcmc_split_lengths.sum(axis=1)

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

    # Note: in the following arrays, the particles are laid out along axis 0 and the
    # splits are laid out along axis 1.
    # Now we need to set things up to translate between split indexing and branch
    # indexing.
    # The ith entry of this array gives the index of the split corresponding to the ith
    # branch.
    branch_to_split = np.array(inst.get_psp_indexer_representations()[0][0])
    # The ith entry of this array gives the index of the branch corresponding to the ith
    # split.
    split_to_branch = np.empty_like(branch_to_split)
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

    def log_like_or_grad_with(split_lengths, grad=False):
        """Calculate log likelihood or the gradient with given split
        lengths."""
        branch_lengths[:] = translate_splits_to_branches(split_lengths)
        if grad:
            _, log_grad = inst.branch_gradients()[0]
            # This :-2 is because of the two trailing zeroes that appear at the end of
            # the gradient.
            result = translate_branches_to_splits(np.array(log_grad)[:-2])
        else:
            result = np.array(inst.log_likelihoods())[0]
        return result

    def phylo_log_like(split_lengths_arr):
        """Calculate phylogenetic log likelihood for each of the split length
        assignments laid out along axis 1."""
        return np.apply_along_axis(log_like_or_grad_with, 1, split_lengths_arr)

    def grad_phylo_log_like(split_lengths_arr):
        return np.apply_along_axis(
            lambda x: log_like_or_grad_with(x, grad=True), 1, split_lengths_arr
        )

    def log_exp_prior(x, rate=10):
        return np.log(rate) - np.sum(rate * x, axis=1)

    def grad_log_exp_prior(x, rate=10):
        return -rate

    def phylo_log_upost(split_lengths_arr):
        """The unnormalized phylogenetic posterior with an Exp(10) prior."""
        return phylo_log_like(split_lengths_arr) + log_exp_prior(split_lengths_arr)

    def grad_phylo_log_upost(split_lengths_arr):
        """The unnormalized phylogenetic posterior with an Exp(10) prior."""
        return grad_phylo_log_like(split_lengths_arr) + grad_log_exp_prior(
            split_lengths_arr
        )

    model = vip.continuous_models.of_name(
        model_name, variable_count=len(branch_lengths), particle_count=particle_count
    )
    model.mode_match(last_sampled_split_lengths)
    opt = vip.optimizers.of_name("bump", model)
    start_time = timeit.default_timer()
    opt.gradient_steps(phylo_log_upost, grad_phylo_log_upost, step_count)
    gradient_time = timeit.default_timer() - start_time
    opt_trace = pd.DataFrame({"elbo": opt.trace}).reset_index()

    fit_sample = pd.DataFrame(model.sample(len(mcmc_split_lengths)))
    fit_sample["total"] = fit_sample.sum(axis=1)
    fit_sample["type"] = "vb"
    mcmc_split_lengths["type"] = "mcmc"
    fitting_results = pd.concat(
        [fit_sample.melt(id_vars="type"), mcmc_split_lengths.melt(id_vars="type")]
    )
    fitting_results["variable"] = fitting_results["variable"].astype(str)
    final_elbo = model.elbo_estimate(
        phylo_log_upost, particle_count=particle_count_for_final_elbo_estimate
    )

    run_details = {"gradient_time": gradient_time, "final_elbo": final_elbo}

    return run_details, opt_trace, fitting_results
