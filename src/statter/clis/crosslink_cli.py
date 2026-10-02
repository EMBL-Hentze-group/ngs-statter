import os

import rich_click as click

from statter.analyzers.roi_xlink_wrapper import RoiOverlapWrapper
from statter.clis.cli_utils import (
    fig_options,
    normalizer,
    overlap_options,
    run_options,
)
from statter.parsers.csv_meta_parser import MetaReader
from statter.plotters.roi_xlink_heatmap import HeatMapper
from statter.plotters.roi_xlink_line_plotter import LinePlotter

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])  # noqa: C408


@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.version_option()
def crosslink() -> None:
    """
    Ops for overlapping crosslinking data over secondary structure/primary motif regions
    """


@crosslink.command("crosslink-meta-example")
def crosslink_meta_example() -> None:
    """
    Print example CSV metadata file for crosslink plotting
    """
    MetaReader.xlink_example()


@crosslink.command(
    "count-crosslinks", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@overlap_options
@click.option(
    "--norm",
    "norm",
    default="cpm",
    help="Normalization method: 'none' or 'cpm'[Counts per million]",
    type=click.Choice(["none", "cpm"]),
    show_default=True,
    callback=normalizer,
)
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
    xov = RoiOverlapWrapper(
        meta_data=csv,
        region=bed,
        l=l,
        r=r,
        unstranded=unstranded,
        most_5prime=most_5prime,
        smoothing_window=smoothing_window,
        out_file=out_table,
        tmp_dir=tmpdir,
    )
    xov.crosslink_overlap(norm_method=norm)


@crosslink.command(
    "crosslink-line-plot", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@overlap_options
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
    xov = RoiOverlapWrapper(
        meta_data=csv,
        region=bed,
        l=l,
        r=r,
        unstranded=unstranded,
        most_5prime=most_5prime,
        smoothing_window=smoothing_window,
        out_file=out_table,
        tmp_dir=tmpdir,
    )
    xov.crosslink_overlap(norm_method=norm)
    region_max_len = xov.region_max_len
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
@overlap_options
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
    xov = RoiOverlapWrapper(
        meta_data=csv,
        region=bed,
        l=l,
        r=r,
        unstranded=unstranded,
        most_5prime=most_5prime,
        smoothing_window=smoothing_window,
        out_file=out_table,
        tmp_dir=tmpdir,
    )
    xov.crosslink_overlap(norm_method=norm)
    region_max_len = xov.region_max_len
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
