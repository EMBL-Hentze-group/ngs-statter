import re
from dis import show_code
from email.policy import default
from multiprocessing import context
from typing import Optional

import click
from numpy import require

from statter.parsers.stats_collector import (
    SampleStats,
    compile_all_sample_stats,
)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--sample_name",
    "-n",
    "sample_name",
    required=True,
    help="Sample name",
    type=str,
)
@click.option(
    "--raw_reads",
    "-r",
    "raw_reads",
    required=True,
    help="Path to seqkit stats output for raw reads",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--first_trim",
    "-f",
    "first_trim",
    required=True,
    help="Path to seqkit stats output for first trimming step",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--second_trim",
    "-s",
    "second_trim",
    help="Path to seqkit stats output for second trimming step [Optional]",
    required=False,
    default=None,
    show_default=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--rRNA_mapped",
    "-m",
    "rRNA_mapped",
    help="Path to seqkit stats output for rRNA mapped reads [Optional]",
    required=False,
    default=None,
    show_default=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--rRNA_free",
    "-u",
    "rRNA_free",
    help="Path to seqkit stats output for rRNA free reads [Optional]",
    required=False,
    default=None,
    show_default=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--align_stats",
    "-a",
    "align",
    required=True,
    help="Path to bam statistics json file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--dedup_stats",
    "-d",
    "dedup",
    help="Path to bam statistics json file after UMI deduplication [Optional]",
    required=False,
    default=None,
    show_default=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--output",
    "-o",
    "output",
    required=True,
    help="Output json file for collected statistics",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
)
def all_stats(
    sample_name,
    raw_reads,
    first_trim,
    second_trim,
    rRNA_mapped,
    rRNA_free,
    align,
    dedup,
    output,
):
    """all_stats
    Collect statistics from various steps and compile into a single json file

    Args:
        sample_name: Sample name
        raw_reads: Path to seqkit stats output for raw reads
        first_trim: Path to seqkit stats output for first trimming step
        second_trim: [Optional] Path to seqkit stats output for second trimming step
        rRNA_mapped: [Optional] Path to seqkit stats output for rRNA mapped reads
        rRNA_free: [Optional] Path to seqkit stats output for rRNA free reads
        align: Path to bam statistics json file
        dedup: [Optional] Path to bam statistics json file after UMI deduplication
        output: Output json file for collected statistics
    """
    stats_collector = SampleStats(
        sample=sample_name,
        raw=raw_reads,
        first_trim=first_trim,
        align=align,
        out=output,
        second_trim=second_trim,
        rRNA_free=rRNA_free,
        rRNA_mapped=rRNA_mapped,
        dedup=dedup,
    )
    stats_collector.to_json()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--output",
    "-o",
    "output",
    required=True,
    help="Output file name for all sample statistics [tsv format]",
    type=click.Path(file_okay=True, dir_okay=False, writable=True),
)
@click.argument(
    "files",
    nargs=-1,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    # help="List of files to process. [output json files from `all_stats -h`]",
)
def compile_stats(output, files):
    """
    Collect statistics from multiple samples and compile into a single json file

    Args:
        files: List of files to process. [output json files from `all_stats -h`]
    """
    compile_all_sample_stats(output, list(files))
