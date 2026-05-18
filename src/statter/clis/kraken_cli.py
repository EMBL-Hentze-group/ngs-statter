import rich_click as click

from statter.parsers.kraken_report_parser import Kraken2
from statter.parsers.taxonomy_parser import Taxonomy

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def kraken() -> None:
    """
    Kraken2 processing utilities
    """


@kraken.command(
    "collect-reports", context_settings=CONTEXT_SETTINGS, no_args_is_help=True
)
@click.argument(
    "reports",
    required=False,
    # help="List of Kraken2 report files (see kraken2 -h). Either --kraken_dir or --reports MUST be provided",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    nargs=-1,
)
@click.option(
    "--nodes",
    "nodes",
    required=True,
    help="NCBI Taxonomy db nodes.dmp file, see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--names",
    "names",
    required=True,
    help="NCBI Taxonomy db names.dmp file, see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
@click.option(
    "--kraken-dir",
    "kraken_dir",
    required=False,
    help="Directory containing Kraken2 classification reports (see https://github.com/DerrickWood/kraken2). Either --kraken-dir or reports as space separated arguments MUST be provided",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option(
    "--pattern",
    "pattern",
    default="*.txt",
    help="Kraken report naming pattern. Used only if --kraken-dir is provided",
    show_default=True,
    type=str,
)
@click.option(
    "--out-csv",
    "output_csv",
    required=True,
    help="Output csv file name",
    type=click.Path(exists=False),
)
def kraken_report_aggregator(
    nodes: str,
    names: str,
    kraken_dir: str,
    pattern: str,
    reports: list[str],
    output_csv: str,
):
    """
    Collect kraken2 reports into a single file\b

    Aggregate Kraken2 species classification report. Given NCBI taxonomy database nodes.dmp file names.dmp file, either a directory with Kraken2 classification results
    or a list of Kraken2 report files, aggregate all classification results into one file. Either --kraken-dir option or kraken reports as space separated arguments must be provided.
    """
    if kraken_dir is None and len(reports) == 0:
        raise RuntimeError(
            "Either --kraken-dir or kraken reports as space separated arguments must be provided"
        )
    if kraken_dir is not None and len(reports) > 0:
        raise RuntimeError(
            "Only one of --kraken-dir or kraken reports as space separated arguments should be provided"
        )
    tax = Taxonomy(nodes=nodes, names=names)
    taxonomy = tax.build_dag()
    k2_agg = Kraken2(
        taxonomy=taxonomy,
        output_file=output_csv,
        report_folder=kraken_dir,
        report_files=reports,
        pattern=pattern,
    )
    k2_agg.aggregate_reports()
