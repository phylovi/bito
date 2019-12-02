"""
Methods for running and benchmarking performance.

Documentation in the command line interface.
"""

import os
import timeit

import libsbn
import numpy as np
import pandas as pd

import vip.burrito


# Documentation is in the CLI.
# `*` forces everything after to be keyword-only.
def fixed(
    data_path,
    *,
    branch_model_name,
    scalar_model_name,
    optimizer_name,
    step_count,
    particle_count,
    thread_count
):
    data_path = os.path.normpath(data_path)
    data_id = os.path.basename(data_path)
    mcmc_nexus_path = os.path.join(data_path, data_id + "_out.t")
    fasta_path = os.path.join(data_path, data_id + ".fasta")
    burn_in_fraction = 0.1
    particle_count_for_final_elbo_estimate = 10000
    phylo_model_specification = libsbn.PhyloModelSpecification(
        substitution="JC69", site="constant", clock="strict"
    )
    # Read MCMC run and get split lengths.
    mcmc_inst = libsbn.instance("mcmc_inst")
    mcmc_inst.read_nexus_file(mcmc_nexus_path)
    burn_in_count = int(burn_in_fraction * mcmc_inst.tree_count())
    mcmc_inst.tree_collection.erase(0, burn_in_count)
    mcmc_inst.process_loaded_trees()
    ragged = [np.array(a) for a in mcmc_inst.split_lengths()]
    mcmc_split_lengths = pd.concat(
        [pd.DataFrame({"variable": idx, "value": a}) for idx, a in enumerate(ragged)],
        sort=False,
    )
    last_sampled_split_lengths = np.array([a[-1] for a in ragged])

    burro = vip.burrito.Burrito(
        mcmc_nexus_path=mcmc_nexus_path,
        burn_in_fraction=burn_in_fraction,
        fasta_path=fasta_path,
        phylo_model_specification=phylo_model_specification,
        branch_model_name=branch_model_name,
        scalar_model_name=scalar_model_name,
        optimizer_name=optimizer_name,
        particle_count=particle_count,
        thread_count=thread_count,
    )
    burro.branch_model.mode_match(last_sampled_split_lengths)

    start_time = timeit.default_timer()
    burro.gradient_steps(step_count)
    gradient_time = timeit.default_timer() - start_time
    opt_trace = pd.DataFrame({"elbo": burro.opt.trace}).reset_index()

    # We sample from our fit model as many times as there were trees in our MCMC sample.
    fit_sample = pd.DataFrame(burro.branch_model.sample_all(mcmc_inst.tree_count()))
    fit_sample["type"] = "vb"
    mcmc_split_lengths["type"] = "mcmc"
    fitting_results = pd.concat(
        [fit_sample.melt(id_vars="type"), mcmc_split_lengths], sort=False
    )
    fitting_results["variable"] = fitting_results["variable"].astype(str)
    final_elbo = burro.estimate_elbo(
        particle_count=particle_count_for_final_elbo_estimate
    )

    run_details = {"gradient_time": gradient_time, "final_elbo": final_elbo}

    return run_details, opt_trace, fitting_results
