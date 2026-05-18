import rich_click as click

from statter.parsers.bam_parser import BamParser

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.version_option()
def alignment() -> None:
    """
    Collect run statistics from alignments
    """


@alignment.command("bam", context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.option(
    "--bam",
    "bam",
    required=True,
    help="Alignments BAM file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--min-q",
    "min_q",
    default=0,
    help="Minimum alignment quality",
    show_default=True,
    type=int,
)
@click.option(
    "--out-json",
    "out_json",
    required=True,
    help="Output file to write json formatted data",
    type=click.Path(exists=False),
)
def basic_alignment_stats(bam: str, min_q: int, out_json: str) -> None:
    """
    Compute basic alignment statistics for a generic BAM file.

    \b
    The basic alignment statistics include the following metrics:
    - Input reads: The total number of reads in the BAM file.
    - Mapped: The number of reads that are mapped to the reference genome.
    - Unmapped: The number of reads that are not mapped to the reference genome.

    The alignment statistics are written to the specified output JSON file.
    """
    bam_parser = BamParser(bam=bam, min_q=min_q, ignore_duplicate=False)
    bam_parser.alignment_stats_rs(out_json=out_json)


@alignment.command("STAR", context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.option(
    "--bam",
    "bam",
    required=True,
    help="Alignments BAM file (MUST be co-ordinate sorted and indexed)",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--min-q",
    "min_q",
    default=0,
    help="Minimum alignment quality",
    show_default=True,
    type=int,
)
@click.option(
    "--out-json",
    "out_json",
    required=True,
    help="Output file to write json formatted data",
    type=click.Path(exists=False),
)
def alignment_stats_STAR(bam: str, min_q: int, out_json: str) -> None:
    """
    Compute alignment statistics for a STAR aligned BAM file and output as JSON.
    """
    bam_parser = BamParser(bam=bam, min_q=min_q, ignore_duplicate=False)
    bam_parser.STAR_alignment_stats_rs(out_json=out_json)
