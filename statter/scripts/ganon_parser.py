import warnings
from logging import warning

import click

from statter.parsers.ganon_abundance_parser import GanonTreParser
from statter.parsers.ganon_unclassified_fastq_writer import GanonUnclassifiedFq
from statter.plotters.ganon_contamination_vs_unclassified_read_length_distribution import (
    GanonReadLengthPlot,
)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

"""
Write ganon unclassified reads into a separate fastq file and 
write unclassified reads read length distribution to json
"""


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--unc",
    "unc",
    required=True,
    help="Ganon unclassifed read ids output file, see 'ganon classify -h' for details",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--src",
    "src_fastq",
    required=True,
    help="Source fastq file, used as input to 'ganon classify'",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--unclassified",
    "target_fastq",
    required=True,
    help="File name to write ganon unclassified reads in fastq format (supports .gz files)",
)
@click.option(
    "--read_length",
    "rl_json",
    required=True,
    help="Write unclassified reads read length distribuion to json",
)
def unclassifed_read_parser(
    unc: str, src_fastq: str, target_fastq: str, rl_json
) -> None:
    """
    Write ganon2 unclassified reads into a separate fastq file and compute read length distribution\b

    **Warning**: This script is deprecated
    """
    warnings.warn(
        "This script is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    ganon_unc = GanonUnclassifiedFq(
        unc=unc, src_fq=src_fastq, target_fq=target_fastq, read_length_json=rl_json
    )
    ganon_unc.write_unclassified_reads()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--json_dir",
    "json_dir",
    required=True,
    help="Directory containing json formatted read_length files after ganon classification (see ganon_unc_to_fastq -h)",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
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
    "--html",
    "output_html",
    required=True,
    help="Plot output file name",
    type=click.Path(exists=False),
)
def ganon_read_length_plotter(
    json_dir: str, output_html: str, pattern: str = "*.json", nsamples: int = 1
) -> None:
    """
    Plot read lengths after ganon classification using json formatted read length distribution files\b

    **Warning**: This script is deprecated
    """
    warnings.warn(
        "This script is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    grlp = GanonReadLengthPlot(
        json_folder=json_dir,
        output_file=output_html,
        pattern=pattern,
        nsamples=nsamples,
    )
    grlp.plot()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--result_dir",
    "ganon_result_dir",
    required=True,
    help="Directory containing ganon species abundance estimation files (*.tre) after ganon classification",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
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
    "--out",
    "output_csv",
    required=True,
    help="Output csv file name",
    type=click.Path(exists=False),
)
def ganon_abundance_aggregator(
    ganon_result_dir: str,
    output_csv: str,
    pattern: str = "*.tre",
) -> None:
    """
    collect ganon report (*.tre) files and aggregate results across multiple samples\b

    **Warning**: This script is deprecated
    """
    warnings.warn(
        "This script is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    gtp = GanonTreParser(
        result_folder=ganon_result_dir,
        output_file=output_csv,
        pattern=pattern,
    )
    gtp.aggregate_results()
