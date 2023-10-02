from parsers.fastq_stats import FqLength
from plotters.read_length_distribution import ReadLengthPlot
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
    For a given fastq file, compute read lenghts and write the distribution to json format
    """
    flp = FqLength(fq=fastq, out_file=json)
    flp.read_length()


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
@click.option("--html", "output_html", required=True, help="Plot output file name")
def read_length_plotter(
    json_dir: str, output_html: str, pattern: str, nsamples: int = 1
) -> None:
    """
    Read in  json formatted files and output html plot
    """
    rlp = ReadLengthPlot(
        json_folder=json_dir,
        output_file=output_html,
        min_read_length=min_read_length,
        pattern=pattern,
        nsamples=nsamples,
    )
    rlp.plot()
