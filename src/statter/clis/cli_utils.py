import os
from functools import wraps
from warnings import warn

import rich_click as click


def to_inches(ctx, param, value):
    """
    Convert centimeters to inches for figure dimensions, as matplotlib uses inches for figure size. This callback function is used for the --fig-width and --fig-height options.
    """
    if value is None:
        return None
    try:
        return value / 2.54
    except TypeError:
        raise click.BadParameter("Value must be a number")


def set_threads(ctx, param, value):
    """
    Validate the number of threads requested by user
    """
    ncpus: int = os.cpu_count()  # type: ignore
    if value > ncpus:
        warn(
            f"Requested {value} threads, but only {ncpus} CPUs available. Using {max(1, ncpus - 1)} threads.",
            RuntimeWarning,
        )
        value = max(1, ncpus - 1)
    return value


def normalizer(ctx, param, value):
    """
    Normalize the input value. If the value is "none", return None.
    """
    if value == "none":
        return None
    return value


def all_to_none(ctx, param, value):
    """
    Convert the string "all" to None. Useful for options where "all" means no filtering.
    """
    if value == "all":
        return None
    return value


def vcf_norm_options(func):
    @click.option(
        "--genome-fasta",
        "genome_fasta",
        required=False,
        default=None,
        help="Path to the genome FASTA file",
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
        show_default=True,
    )
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def vcf_query_options(func):
    @click.option(
        "--variant-type",
        "variant_type",
        default="all",
        help="Type of variant to query (e.g., snp, indel,...)",
        type=click.Choice(["all", "snp", "mnp", "indel", "insertion", "deletion"]),
        show_default=True,
        callback=all_to_none,
    )
    @click.option(
        "--read-depth",
        "read_depth",
        required=False,
        default=0,
        help="Minimum read depth for querying variants",
        type=click.IntRange(min=0),
        show_default=True,
    )
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def overlap_options(func):
    @click.option(
        "--metadata",
        "csv",
        required=True,
        help="CSV metadata file specifying crosslinking site files and sample information (see `ngs-statter crosslink-meta-example` OR `ngs-statter vcf-meta-example`)",
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
    )
    @click.option(
        "--bed",
        "bed",
        required=True,
        help="BED file specifying secondary structure/ primary motif regions of interests (supports .gz files)",
        type=click.Path(exists=True, file_okay=True, dir_okay=False),
    )
    @click.option(
        "--out-table",
        "out_table",
        required=True,
        help="Output file to write the aggregated table (always .parquet format)",
        type=click.Path(exists=False),
    )
    @click.option(
        "--l",
        "l",
        default=100,
        help="5' extension length for regions in BED file",
        show_default=True,
        type=click.IntRange(min=0),
    )
    @click.option(
        "--r",
        "r",
        default=100,
        help="3' extension length for regions in BED file",
        show_default=True,
        type=click.IntRange(min=0),
    )
    @click.option(
        "--unstranded",
        "unstranded",
        is_flag=True,
        show_default=True,
        default=False,
        help="If this flag is set, ignore strand information in the BED file and treat all regions as unstranded",
    )
    @click.option(
        "--most-5prime",
        "most_5prime",
        is_flag=True,
        show_default=True,
        default=False,
        help="If bed regions overlap, only keep the most 5' region out of the overlapping regions",
    )
    @click.option(
        "--sw",
        "smoothing_window",
        help="When plotting smooth crosslink sites using moving average. Use these many adjacent bases to compute moving average",
        type=click.IntRange(min=1),
        default=5,
        show_default=True,
    )
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def run_options(func):
    @click.option(
        "--tmpdir",
        "tmpdir",
        default=None,
        help="Temporary directory to use (default: system temp folder)",
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
        show_default=True,
    )
    @click.option(
        "--threads",
        "threads",
        default=4,
        help="Number of threads to use",
        type=click.IntRange(min=1),
        callback=set_threads,
        show_default=True,
    )
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def fig_options(func):
    @click.option(
        "--fig-width",
        "width",
        default=30,
        help="Figure width in centimeters",
        type=click.FloatRange(min=1),
        callback=to_inches,
        show_default=True,
    )
    @click.option(
        "--fig-height",
        "height",
        default=27,
        help="Figure height in centimeters",
        type=click.FloatRange(min=1),
        callback=to_inches,
        show_default=True,
    )
    @click.option(
        "--xlabel",
        "xlabel",
        default="X axis label",
        help="X axis label for the plot",
        type=str,
        show_default=True,
    )
    @click.option(
        "--ylabel",
        "ylabel",
        default="Y axis label",
        help="Y axis label for the plot",
        type=str,
        show_default=True,
    )
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper
