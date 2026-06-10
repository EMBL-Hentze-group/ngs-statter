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


@click.group(
    name="ngs-statter",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
)
@click.version_option(version=version("ngs-statter"), prog_name="ngs-statter")
def runner() -> None:
    """
    for command specific help, use: **ngs-statter *command* -h**

    The list of available commands are shown below grouped by their source.
    """
    pass


sources = [alignment, crosslink, fastq, flexbar, genetype, kraken, sample]

click.rich_click.COMMAND_GROUPS["ngs-statter"] = []
for source in sources:
    for name, cmd in source.commands.items():
        runner.add_command(cmd, name=name)  # add the command to the main cli group
    group_panel = {
        "name": f"{source.name} commands",
        "commands": list(source.commands.keys()),
    }
    click.rich_click.COMMAND_GROUPS["ngs-statter"].append(group_panel)
