import os

import rich_click as click

from statter.analyzers.roi_xlink_wrapper import RoiOverlapWrapper
from statter.clis.cli_utils import (
    fig_options,
    overlap_options,
    run_options,
    vcf_norm_options,
    vcf_query_options,
)
from statter.parsers.csv_meta_parser import MetaReader
from statter.plotters.roi_xlink_line_plotter import LinePlotter

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])  # noqa: C408


@click.group(context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
@click.version_option()
def vcf() -> None:
    """
    Ops for overlapping CIMS (Crosslink Induced Mutation Sites) data over secondary structure/primary motif regions
    """


@vcf.command("vcf-meta-example")
def vcf_meta_example() -> None:
    """
    Print example CSV metadata file for VCF plotting
    """
    MetaReader.vcf_example()


@vcf.command(
    "parse-alt-allele-depth", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@overlap_options
@vcf_query_options
@vcf_norm_options
@run_options
def parse_alt_allele_depth(
    csv,
    genome_fasta,
    variant_type,
    read_depth,
    bed,
    out_table,
    l,
    r,
    unstranded,
    most_5prime,
    smoothing_window,
    tmpdir,
    threads,
) -> None:
    """
    Parse alternate allele depth from VCF files based on the provided options and save to a table.
    """
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    # Implementation goes here
    ad_ov = RoiOverlapWrapper(
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
    ad_ov.allele_depth_overlap(
        genome_fa=genome_fasta,
        variant_type=variant_type,
        dp=read_depth,
        ncores=threads,
    )


@vcf.command(
    "plot-alt-allele-depth", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@overlap_options
@vcf_query_options
@vcf_norm_options
@click.option(
    "--out-fig",
    "out_fig",
    required=True,
    help="Output file to write the plot (svg/pdf/png)",
    type=click.Path(exists=False),
)
@fig_options
@click.option(
    "--title",
    "title",
    default="Crosslink Induced Mutations: Alternate allele depth",
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
def plot_alt_allele_depth(
    csv,
    genome_fasta,
    variant_type,
    read_depth,
    bed,
    out_table,
    l,
    r,
    unstranded,
    most_5prime,
    smoothing_window,
    out_fig,
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
    Parse alternate allele depth from VCF files based on the provided options and plot allele depths as line plots
    """
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    ad_ov = RoiOverlapWrapper(
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
    ad_ov.allele_depth_overlap(
        genome_fa=genome_fasta,
        variant_type=variant_type,
        dp=read_depth,
        ncores=threads,
    )
    region_max_len = ad_ov.region_max_len
    lp = LinePlotter(out_table)
    lp.plot(
        output=out_fig,
        width=width,
        height=height,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        ymax=ymax,
        show_group_mean=show_group_mean,
        errorbar=errorbar,
        roi_length=region_max_len,
    )
