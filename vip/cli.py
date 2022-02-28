"""The ``vip`` command line interface."""
import pprint
import click


@click.group()
def cli_benchmark():
    pass


@cli_benchmark.command(name="benchmark")
@click.option(
    "--branch-model",
    type=click.Choice(["split", "psp"]),
    default="split",
    show_default=True,
)
@click.option(
    "--scalar-model",
    type=click.Choice(
        ["lognormal", "tf_lognormal", "tf_truncated_lognormal", "tf_gamma"]
    ),
    default="lognormal",
    show_default=True,
)
@click.option(
    "--optimizer",
    type=click.Choice(["simple", "bump"]),
    default="simple",
    show_default=True,
)
@click.option(
    "--step-count",
    default=5,
    help="Number of gradient descent steps to take.",
    show_default=True,
)
@click.option(
    "--particle-count",
    default=10,
    help="Number of particles to use for stochastic gradient estimation.",
    show_default=True,
)
@click.option(
    "--thread-count", default=4, help="Number of threads to use.", show_default=True
)
@click.option(
    "--out-prefix", default=None, help="Path prefix to which output should be saved."
)
@click.argument("data-path")
def benchmark(
    branch_model,
    scalar_model,
    optimizer,
    step_count,
    particle_count,
    thread_count,
    out_prefix,
    data_path,
):
    """Do a benchmarking comparison to an MCMC run.

    ``DATA_PATH`` is a path to a directory, which say is named ``X``.

    We assume that ``X`` contains:

    ``X_out.t``
      an MCMC run on a fixed tree topology, and

    ``X.fasta``
      a FASTA file with the same sequence data as used for MCMC.
    """
    # This splendid non-Pythonic import in a function means that the CLI is fast unless
    # we do something.
    import vip.benchmark

    print("Starting validation:")
    pprint.pprint(locals())

    run_details, opt_trace, fitting_results = vip.benchmark.fixed(
        data_path,
        branch_model_name=branch_model,
        scalar_model_name=scalar_model,
        optimizer_name=optimizer,
        step_count=step_count,
        particle_count=particle_count,
        thread_count=thread_count,
    )
    if out_prefix is not None:
        opt_trace.to_csv(out_prefix + "_opt_trace.csv")
        fitting_results.to_csv(out_prefix + "_fitting_results.csv")
    pprint.pprint(run_details)


@click.group()
def cli_dag_to_dot():
    pass


@cli_dag_to_dot.command(name="dag-to-dot")
@click.option(
    "-fasta",
    "--fasta-path",
    required=True,
    prompt=True,
    type=click.Path(exists=True),
    help="File path to fasta file.",
)
@click.option(
    "-newick",
    "--newick-path",
    required=True,
    prompt=True,
    type=click.Path(exists=True),
    help="File path to newick file.",
)
@click.option(
    "-output",
    "--output-path",
    required=True,
    prompt=True,
    type=click.Path(),
    help="Path to write output dot and SVG file.",
)
@click.option(
    "-edges",
    "--edge-labels",
    default=False,
    show_default=True,
    help="Specify whether or not edge labels will appear in output.",
)
def dag_to_dot(fasta_path, newick_path, output_path, edge_labels):
    """Convert a subsplit DAG to a .dot and .svg file.

    Provide file paths relative to the directory where the command is run.

    The command will output the .dot file and render it with graphviz into a .svg file.
    """
    import bito
    import tempfile
    import graphviz

    mmap_file = tempfile.mkstemp(suffix=".data")[1]

    inst = bito.gp_instance(mmap_file)
    inst.read_fasta_file(fasta_path)
    inst.read_newick_file(newick_path)
    inst.make_engine()
    inst.subsplit_dag_to_dot(output_path, edge_labels)

    graphviz.render("dot", "svg", output_path)


cli = click.CommandCollection(sources=[cli_benchmark, cli_dag_to_dot])
