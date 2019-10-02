import click

@click.group()
def cli():
    pass

@cli.command()
def benchmark():
    # This splendid non-Pythonic import in a function means that the CLI is fast unless
    # we do something.
    import vip.benchmark
    vip.benchmark.fixed()
