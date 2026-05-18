import os
import warnings
from functools import wraps

import rich_click as click

from statter.analyzers.roi_xlink_overlapper import RegionXlinkOverlapFinder
from statter.parsers.csv_meta_parser import MetaReader
from statter.plotters.roi_xlink_heatmap import HeatMapper
from statter.plotters.roi_xlink_line_plotter import LinePlotter

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


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
        warnings.warn(
            f"Requested {value} threads, but only {ncpus} CPUs available. Using {max(1, ncpus - 1)} threads."
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


def crosslink_options(func):
    @click.option(
        "--metadata",
        "csv",
        required=True,
        help="CSV metadata file specifying crosslinking site files and sample information (see `statter csv-meta-example`)",
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
        "--norm",
        "norm",
        default="cpm",
        help="Normalization method: 'none' or 'cpm'[Counts per million]",
        type=click.Choice(["none", "cpm"]),
        show_default=True,
        callback=normalizer,
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
        default="Relative crosslink positions",
        help="X axis label for the plot",
        type=str,
        show_default=True,
    )
    @click.option(
        "--ylabel",
        "ylabel",
        default="Crosslink counts",
        help="Y axis label for the plot",
        type=str,
        show_default=True,
    )
    @wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.version_option()
def crosslink() -> None:
    """
    Plot crosslinking data over secondary structure/primary motif regions
    """


@crosslink.command("csv-meta-example")
def csv_meta_example() -> None:
    """
    Print example CSV metadata file for crosslink plotting
    """
    MetaReader.metadata_example()


@crosslink.command(
    "count-crosslinks", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@crosslink_options
@run_options
def count_crosslinks(
    csv,
    bed,
    out_table,
    l,
    r,
    most_5prime,
    norm,
    unstranded,
    smoothing_window,
    tmpdir,
    threads,
) -> None:
    """
    Count crosslinking sites over secondary structure/primary motif regions.
    """
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    with RegionXlinkOverlapFinder(
        metadata=csv,
        region=bed,
        l=l,
        r=r,
        unstranded=unstranded,
        most_5prime=most_5prime,
        norm_method=norm,
        smoothing_window=smoothing_window,
        tmpdir=tmpdir,
    ) as overlapper:
        overlapper.find_overlaps(out_path=out_table)


@crosslink.command(
    "crosslink-line-plot", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@crosslink_options
@click.option(
    "--out-fig",
    "out",
    required=True,
    help="Output file to write the plot (svg/pdf/png)",
    type=click.Path(exists=False),
)
@fig_options
@click.option(
    "--title",
    "title",
    default="Crosslink profile",
    help="Title for the plot",
    type=str,
    show_default=True,
)
@click.option(
    "--ymax",
    "ymax",
    default=None,
    help="Maximum value for crosslink counts on y axis (determined from data if not set)",
    type=click.FloatRange(min=0.0),
    show_default=True,
)
@click.option(
    "--show-group-mean",
    "show_group_mean",
    is_flag=True,
    default=False,
    help="Whether to show group mean in the plot",
    show_default=True,
)
@click.option(
    "--errorbar",
    "errorbar",
    default=None,
    help="Error bar to use",
    type=click.Choice(["sd", "ci", "pi", "se", "sd", "None"]),
    show_default=True,
)
@run_options
def crosslink_line_plot(
    csv,
    bed,
    out,
    out_table,
    l,
    r,
    most_5prime,
    norm,
    unstranded,
    smoothing_window,
    xlabel,
    ylabel,
    title,
    ymax,
    width,
    height,
    show_group_mean,
    errorbar,
    tmpdir,
    threads,
) -> None:
    """
    Plot crosslinking sites over secondary structure/primary motif regions as line plots.
    """
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    with RegionXlinkOverlapFinder(
        metadata=csv,
        region=bed,
        l=l,
        r=r,
        unstranded=unstranded,
        most_5prime=most_5prime,
        norm_method=norm,
        smoothing_window=smoothing_window,
        tmpdir=tmpdir,
    ) as overlapper:
        overlapper.find_overlaps(out_path=out_table)
        region_max_len = overlapper.region_max_len
    lp = LinePlotter(out_table)
    lp.plot(
        output=out,
        width=width,
        height=height,
        ymax=ymax,
        errorbar=errorbar,
        roi_length=region_max_len,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        show_group_mean=show_group_mean,
    )


@crosslink.command(
    "crosslink-heatmap", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@crosslink_options
@click.option(
    "--out-dir",
    "out",
    required=True,
    help="Output directory for plots. Group specific heatmaps will be written to this directory (svg format)",
    type=click.Path(exists=False),
)
@fig_options
@click.option(
    "--vmin",
    "vmin",
    default=None,
    help="Minimum value for crosslink counts on y axis (determined from data if not set)",
    type=click.FloatRange(min=0.0),
    show_default=True,
)
@click.option(
    "--vmax",
    "vmax",
    default=None,
    help="Maximum value for crosslink counts on y axis (determined from data if not set)",
    type=click.FloatRange(min=1.0),
    show_default=True,
)
@run_options
def crosslink_heatmap(
    csv,
    bed,
    out,
    out_table,
    l,
    r,
    most_5prime,
    norm,
    unstranded,
    smoothing_window,
    xlabel,
    ylabel,
    vmin,
    vmax,
    width,
    height,
    tmpdir,
    threads,
) -> None:
    """
    Plot crosslinking sites over secondary structure/primary motif regions as heatmaps.
    """
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    with RegionXlinkOverlapFinder(
        metadata=csv,
        region=bed,
        l=l,
        r=r,
        unstranded=unstranded,
        most_5prime=most_5prime,
        norm_method=norm,
        smoothing_window=smoothing_window,
        tmpdir=tmpdir,
    ) as overlapper:
        overlapper.find_overlaps(out_path=out_table)
        region_max_len = overlapper.region_max_len
    hm = HeatMapper(out_table)
    hm.plot(
        outdir=out,
        width=width,
        height=height,
        vmin=vmin,
        vmax=vmax,
        xlabel=xlabel,
        ylabel=ylabel,
        roi_length=region_max_len,
    )
