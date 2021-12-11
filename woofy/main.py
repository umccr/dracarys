import click
from woofy import dog

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version.__version__)
def cli():
    """Somatic Workflow Postprocessing"""


cli.add_command(dog.woofy)
