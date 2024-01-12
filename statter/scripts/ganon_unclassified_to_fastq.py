import click

from statter.parsers.ganon_unclassified_fastq_writer import GanonUnclassifiedFq

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
    "--unclassifed",
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
    Write ganon2 unclassified reads into a separate fastq file and compute read length distribution
    """
    ganon_unc = GanonUnclassifiedFq(
        unc=unc, src_fq=src_fastq, target_fq=target_fastq, read_length_json=rl_json
    )
    ganon_unc.write_unclassified_reads()
