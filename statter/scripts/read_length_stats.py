from statter.parsers.fastq_stats import FqLength
from statter.plotters.read_length_distribution import ReadLengthPlot
from statter.parsers.fastp_json_parser import Fastp
import click

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

"""
Parse length of the reads from a given fastq file and save to json format
"""


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--fq", "fastq", required=True, help="Fastq file (supports .gz files)")
@click.option(
    "--json", "json", required=True, help="File name for json formatted output file"
)
def read_length_parser(fastq: str, json: str) -> None:
    """
    For a given fastq file, compute read lengths and write the distribution to json format
    """
    flp = FqLength(fq=fastq, out_file=json)
    flp.read_length_rs()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--json_dir",
    "json_dir",
    required=True,
    help="Directory containing json formatted read_length files (see read_length_parser -h)",
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
def read_length_plotter(
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
    rlp = ReadLengthPlot(
        json_folder=json_dir,
        output_file=output_html,
        min_read_length=min_read_length,
        pattern=pattern,
        title=title,
        nsamples=nsamples,
    )
    rlp.plot()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--json_dir",
    "trim_dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing json formatted trimming stats from fastp (see fastp -h)",
)
@click.option(
    "--out_csv",
    "out_csv",
    required=True,
    type=click.Path(exists=False),
    help="File name for aggregated output file",
)
@click.option(
    "--first_trim_pattern",
    "first_trim_pattern",
    default="*first_trim.json",
    help="File naming pattern for stats from first trimming",
    show_default=True,
    type=str,
)
@click.option(
    "--second_trim_pattern",
    "second_trim_pattern",
    default="*second_trim.json",
    help="File naming pattern for stats from second trimming",
    show_default=True,
    type=str,
)
def fastp_stats_aggregator(
    trim_dir: str, out_csv: str, first_trim_pattern: str, second_trim_pattern: str
) -> None:
    """
    collect trimming stats from fastp json formatted files for data after first and second trim and aggregate trimming stats into one single file
    """
    fp = Fastp(
        trim_dir=trim_dir,
        out_csv=out_csv,
        first_trim_pattern=first_trim_pattern,
        second_trim_pattern=second_trim_pattern,
    )
    fp.collect_trimming_stats()
