import re

import click
from numpy import require

from statter.parsers.fastq_stats import FqLength
from statter.plotters.read_length_distribution import ReadLengthPlot

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def fastq() -> None:
    """
    Fastq file processing utilities
    """


@fastq.command("parse_read_length", context_settings=CONTEXT_SETTINGS)
@click.option("--fq", "fastq", required=True, help="Fastq file (supports .gz files)")
@click.option(
    "--json", "json", required=True, help="File name for json formatted output file"
)
def parse_read_length(fastq: str, json: str) -> None:
    """
    Compute fastq read length distribution

    \b

    For a given fastq file, compute read lengths and write the distribution to json format
    """
    flp = FqLength(fq=fastq, out_file=json)
    flp.read_length()


@fastq.command("plot_read_length", context_settings=CONTEXT_SETTINGS)
@click.argument(
    "json",
    nargs=-1,
    required=False,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    # help ="List of json formatted read length files (see parse_read_length -h)",
)
@click.option(
    "--json_dir",
    "json_dir",
    required=False,
    help="Directory containing json formatted read_length files (see read_length_parser -h). Either --json_dir or json formatted files as space separated arguments MUST be provided",
)
@click.option(
    "--min_read_length",
    "min_read_length",
    default=15,
    help="Minimum read length that will be kept",
    show_default=True,
    type=int,
)
@click.option(
    "--pattern",
    "pattern",
    default="*.json",
    help="File naming pattern",
    show_default=True,
    type=str,
)
@click.option(
    "--nsamples",
    "nsamples",
    default=1,
    help="Number of replicates per group",
    show_default=True,
    type=int,
)
@click.option(
    "--title",
    "title",
    default="Read length distribution after adapter trimming",
    help="Plot title",
    show_default=True,
    type=str,
)
@click.option("--html", "output_html", required=True, help="Plot output file name")
def plot_read_length(
    json: list[str],
    json_dir: str,
    output_html: str,
    min_read_length: int,
    pattern: str,
    title: str,
    nsamples: int = 1,
) -> None:
    """
    Read in  json formatted files and output html plot
    """
    if json_dir is None and len(json) == 0:
        raise click.UsageError(
            "Either --json_dir option or json formatted files as space separated arguments must be provided"
        )
    if json_dir is not None and len(json) > 0:
        raise click.UsageError(
            "Either --json_dir option or json formatted files as space separated arguments must be provided, not both"
        )
    rlp = ReadLengthPlot(
        json_folder=json_dir,
        json_files=json,
        output_file=output_html,
        min_read_length=min_read_length,
        pattern=pattern,
        title=title,
        nsamples=nsamples,
    )
    rlp.plot()
