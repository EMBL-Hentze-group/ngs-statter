import click

from statter.parsers.flexbar_fastq_parser import HeaderFixer

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def flexbar() -> None:
    """
    Flexbar FASTQ file processing utilities
    """


@flexbar.command("fix_header", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--in_fq",
    "-i",
    "in_fq",
    required=True,
    help="Input FASTQ file with UMI in header comment section, supports .gz files",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--out_fq",
    "-o",
    "out_fq",
    required=True,
    help="Output FASTQ file with fixed UMI position in header, supports .gz files",
    type=click.Path(exists=False, file_okay=True, dir_okay=False),
)
@click.option(
    "--separator",
    "-s",
    "separator",
    default="_",
    help="Separator used to separate UMI from the rest of the header",
    show_default=True,
    type=str,
)
def fix_flexbar_umi_headers(in_fq: str, out_fq: str, separator: str) -> None:
    """
    Fix flexbar UMI headers.
    \b
    Move UMIs from comment section to read identifier in FASTQ headers produced by Flexbar
    """
    fixer = HeaderFixer(in_fq=in_fq, out_fq=out_fq, separator=separator)
    fixer.fix_umi_header()
