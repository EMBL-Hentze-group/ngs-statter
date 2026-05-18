from importlib.metadata import version

import rich_click as click

from statter.clis.bamstat_cli import alignment
from statter.clis.crosslink_cli import crosslink
from statter.clis.flexbar_cli import flexbar
from statter.clis.fq_cli import fastq
from statter.clis.genetype_cli import genetype
from statter.clis.kraken_cli import kraken
from statter.clis.sample_cli import sample

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="statter", context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.version_option(version=version("statter"), prog_name="statter")
def cli() -> None:
    pass


sources = [alignment, crosslink, fastq, flexbar, genetype, kraken, sample]

click.rich_click.COMMAND_GROUPS["statter"] = []
for source in sources:
    for name, cmd in source.commands.items():
        cli.add_command(cmd, name=name)  # add the command to the main cli group
    group_panel = {
        "name": f"{source.name} commands",
        "commands": list(source.commands.keys()),
    }
    click.rich_click.COMMAND_GROUPS["statter"].append(group_panel)


def main() -> None:
    cli()
