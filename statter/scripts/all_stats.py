import click
from typing import Optional
from statter.parsers.sample_stats import SampleStats, _gather_sample_stats

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--first_trim",
    "first_trim",
    required=True,
    help="Fastp .json output from first trimming step",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--second_trim",
    "second_trim",
    required=True,
    help="Fastp .json output from second trimming step",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--rRNA_free",
    "rRNA_free",
    required=True,
    help="rRNA free read length distribution json file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--rRNA_mapped",
    "rRNA_mapped",
    required=True,
    help="rRNA mapped read length distribution json file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--align_stats",
    "align_stats",
    required=True,
    help="Genome alignment statistics json file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--skip_check_rRNA",
    is_flag=True,
    default=False,
    show_default=True,
    help="By defalult, total number of rRNA free reads is checked against the number input reads for genome alignment. If this flag is set, the check is skipped.",
)
@click.option(
    "--dedup_stats",
    "dedup_stats",
    help="Alignment statistics after UMI deduplicateion json file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--output",
    "output",
    required=True,
    help="Output json file",
    type=click.Path(exists=False),
)
def sample_stats(
    first_trim: str,
    second_trim: str,
    rRNA_free: str,
    rRNA_mapped: str,
    align_stats: str,
    output: str,
    skip_check_rRNA: bool = False,
    dedup_stats: Optional[str] = None,
):
    """
    Gather per sample trimming, rRNA removal and genome alignment statistics into one json file
    """
    stat = SampleStats(
        first_trim=first_trim,
        second_trim=second_trim,
        rRNA_free=rRNA_free,
        rRNA_mapped=rRNA_mapped,
        align_stats=align_stats,
        check_rRNA=not skip_check_rRNA,  # default is to check rRNA free reads agaist genome mapped reads
        dedup_stats=dedup_stats,
    )
    stat.collect_stats(output)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--stats_dir",
    "stats_dir",
    required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Directory containing json formatted sample stats stats (see sample_stats -h)",
)
@click.option(
    "--out_csv",
    "out_csv",
    required=True,
    type=click.Path(exists=False),
    help="File name for aggregated output file",
)
@click.option(
    "--suffix",
    "suffix",
    default="*sample_stats.json",
    help="File naming pattern for sample stats",
    show_default=True,
    type=str,
)
def gather_sample_stats(stats_dir: str, out_csv: str, suffix: str) -> None:
    """
    Aggregate all sample stats into one csv file
    \b
    see sample_stats -h
    """
    _gather_sample_stats(stats_dir=stats_dir, out=out_csv, suffix=suffix)
