import click
from dracarys import roar, _version

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(_version.__version__)
def cli():
    """Post-processing DRAGEN workflows."""


cli.add_command(roar.roar)
