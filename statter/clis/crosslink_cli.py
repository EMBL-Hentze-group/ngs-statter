import click

from statter.parsers.crosslink.yaml_parser import YamlCheck
from statter.plotters.crosslink_plotter import Plotter

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def crosslink() -> None:
    """
    Plot crosslinking data over secondary structure/primary motif regions
    """


@crosslink.command("yaml-example")
def yaml_example() -> None:
    """
    Print example YAML file for crosslink plotting
    """
    yc = YamlCheck()
    yc.yaml_example()


@crosslink.command("plot-crosslinks", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--yaml",
    "yaml",
    required=True,
    help="YAML file specifying crosslinking site files and sample information (see `statter yaml-example`)",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--bed",
    "bed",
    required=True,
    help="BED file specifying secondary structure/ primary motif regions",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--out",
    "out",
    required=True,
    help="Output file to write the plot (pdf/png)",
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
    "--norm",
    "norm",
    default="cpm",
    help="Normalization method: 'raw' or 'cpm'[Counts per million]",
    type=click.Choice(["raw", "cpm"]),
    show_default=True,
)
@click.option(
    "--sw",
    "smoothing_window",
    help="When plotting smooth crosslink sites using moving average. Use these many adjacent bases to compute moving average",
    type=click.IntRange(min=1),
    default=5,
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
@click.option(
    "--ymax",
    "ymax",
    default=None,
    help="Maximum value for crosslink counts on y axis(if not set, determined automatically)",
    type=click.IntRange(min=1),
    show_default=True,
)
@click.option(
    "--fig-width",
    "width",
    default=30,
    help="Figure width in centimeters",
    type=click.IntRange(min=10),
    show_default=True,
)
@click.option(
    "--fig-height",
    "height",
    default=27,
    help="Figure height in centimeters",
    type=click.IntRange(min=10),
    show_default=True,
)
@click.option(
    "--errorbar",
    "errorbar",
    default="sd",
    help="Error bar to use",
    type=click.Choice(["sd", "ci", "pi", "se", "sd", "None"]),
    show_default=True,
)
@click.option(
    "--tmpdir",
    "tmpdir",
    default=None,
    help="Temporary directory to use (default: system temp folder)",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    show_default=True,
)
def plot_crosslinks(
    yaml,
    bed,
    out,
    l,
    r,
    norm,
    smoothing_window,
    xlabel,
    ylabel,
    ymax,
    width,
    height,
    errorbar,
    tmpdir,
) -> None:
    """
    Plot crosslinking sites over secondary structure/primary motif regions
    """
    with Plotter(yaml=yaml, region_bed=bed, tmpdir=tmpdir) as pl:
        pl.aggregate_sites(norm=norm)
        pl.index_region(l=l, r=r)
        pl.plot(
            out_file=out,
            smoothing_window=smoothing_window,
            xlabel=xlabel,
            ylabel=ylabel,
            width=width,
            height=height,
            errorbar=errorbar,
            ymax=ymax,
        )
