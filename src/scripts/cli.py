import click

from statter.clis.bamstat_cli import alignment

# from statter.clis.crosslink_cli import crosslink
from statter.clis.flexbar_cli import flexbar
from statter.clis.fq_cli import fastq
from statter.clis.genetype_cli import genetype
from statter.clis.kraken_cli import kraken
from statter.clis.sample_cli import sample

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


clis = click.CommandCollection(
    sources=[flexbar, fastq, alignment, genetype, kraken, sample],
    context_settings=CONTEXT_SETTINGS,
)


def main() -> None:
    clis()


# if groups need to be nested and visible in help output, use the following pattern:
# @click.group(context_settings=CONTEXT_SETTINGS)
# def entry() -> None:
#     """
#     run as: statter <group> <command> [argument] [option]
#     """
#     pass


# entry.add_command(flexbar)
# entry.add_command(fastq)
# entry.add_command(alignment)
# entry.add_command(genetype)
# entry.add_command(kraken)
# entry.add_command(sample)
