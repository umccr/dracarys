import click
import importlib.metadata
from dracarys import roar

__version__ = importlib.metadata.version('dracarys')

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(__version__)
def cli():
    """Post-processing DRAGEN workflows."""


cli.add_command(roar.roar)
