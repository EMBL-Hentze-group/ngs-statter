# from pathlib import Path
import json
import warnings

import click

from statter.statter import single_end_unmapped_fastq

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--bam",
    "bam",
    required=True,
    help="Alignments BAM file (MUST be co-ordinate sorted and indexed)",
)
@click.option(
    "--fq_multi",
    "fq_multi",
    required=True,
    help="File name for multimapped fastq reads (supports .gz files)",
)
@click.option(
    "--fq_other",
    "fq_other",
    required=True,
    help="File name for other unmapped fastq reads (supports .gz files)",
)
@click.option(
    "--stats",
    "stats",
    required=True,
    help="File name for unmapped stats data in json format",
)
def unmapped_fastq_single_end(
    bam: str, fq_multi: str, fq_other: str, stats: str
) -> None:
    """
    parse single end bam file for unmapped reads.\b

    **Warning**: This script is deprecated. Use samtools 1.21 or higher (https://github.com/samtools/samtools/releases/tag/1.21) for the same outputs.\b

    Split unmapped fastqs into two files: (i) all multimappers (ii) other unmapped reads and
    generate unmapped statistics as json
    """
    warnings.warn(
        "Use samtools 1.21 or higher (https://github.com/samtools/samtools/releases/tag/1.21) for the same outputs. This script is deprecated and will be removed in a future release.",
        DeprecationWarning,
        stacklevel=2,
    )
    unmapped_json = {}
    for k, v in single_end_unmapped_fastq(bam, fq_multi, fq_other).items():
        unmapped_json[k] = {"count": v.count, "seq_lens": v.seqlens}
    # if Path(stats).exists()
    with open(stats, "w") as jh:
        json.dump(unmapped_json, jh)
