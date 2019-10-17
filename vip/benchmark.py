import os
import timeit

import libsbn
import numpy as np
import pandas as pd

import vip.burrito


# Documentation is in the CLI.
# `*` forces everything after to be keyword-only.
def fixed(data_path, *, model_name, optimizer_name, step_count, particle_count):
    data_path = os.path.normpath(data_path)
    data_id = os.path.basename(data_path)
    mcmc_nexus_path = os.path.join(data_path, data_id + "_out.t")
    fasta_path = os.path.join(data_path, data_id + ".fasta")
    burn_in_fraction = 0.1
    particle_count_for_final_elbo_estimate = 2000

    # Read MCMC run and get split lengths.
    mcmc_inst = libsbn.instance("mcmc_inst")
    mcmc_inst.read_nexus_file(mcmc_nexus_path)
    mcmc_inst.process_loaded_trees()
    burn_in_count = int(burn_in_fraction * mcmc_inst.tree_count())
    mcmc_split_lengths_np = np.array([np.array(a) for a in mcmc_inst.split_lengths()])
    mcmc_split_lengths = pd.DataFrame(
        mcmc_split_lengths_np[:, burn_in_count:].transpose()
    )
    last_sampled_split_lengths = mcmc_split_lengths.iloc[-1].to_numpy()
    mcmc_split_lengths["total"] = mcmc_split_lengths.sum(axis=1)

    burro = vip.burrito.Burrito(
        mcmc_nexus_path=mcmc_nexus_path,
        fasta_path=fasta_path,
        model_name=model_name,
        optimizer_name=optimizer_name,
        step_count=step_count,
        particle_count=particle_count,
    )
    burro.opt.scalar_model.mode_match(last_sampled_split_lengths)

    start_time = timeit.default_timer()
    burro.gradient_steps(step_count)
    gradient_time = timeit.default_timer() - start_time
    opt_trace = pd.DataFrame({"elbo": burro.opt.trace}).reset_index()

    fit_sample = pd.DataFrame(
        burro.opt.scalar_model.sample(
            len(mcmc_split_lengths),
            which_variables=np.arange(burro.scalar_model.variable_count),
        )
    )
    fit_sample["total"] = fit_sample.sum(axis=1)
    fit_sample["type"] = "vb"
    mcmc_split_lengths["type"] = "mcmc"
    fitting_results = pd.concat(
        [fit_sample.melt(id_vars="type"), mcmc_split_lengths.melt(id_vars="type")]
    )
    fitting_results["variable"] = fitting_results["variable"].astype(str)
    final_elbo = burro.elbo_estimate(
        particle_count=particle_count_for_final_elbo_estimate
    )

    run_details = {"gradient_time": gradient_time, "final_elbo": final_elbo}

    return run_details, opt_trace, fitting_results
