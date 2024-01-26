import click
from statter.parsers.taxonomy_parser import Taxonomy
from statter.parsers.kraken_report_parser import Kraken2

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
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
    "--kraken_dir",
    "kraken_dir",
    required=True,
    help="Directory containing Kraken2 classification reports (see kraken2 -h)",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option(
    "--pattern",
    "pattern",
    default="*.txt",
    help="Kraken report naming pattern",
    show_default=True,
    type=str,
)
@click.option(
    "--out",
    "output_csv",
    required=True,
    help="Output csv file name",
    type=click.Path(exists=False),
)
def kraken_report_aggregator(
    nodes: str, names: str, kraken_dir: str, pattern: str, output_csv: str
):
    """
    Aggregate Kraken2 species classification report
    Given NCBI taxonomy database nodes.dmp file names.dmp file, a directory with Kraken2 classification results
    aggregate all classification results into one file
    """
    tax = Taxonomy(nodes=nodes, names=names)
    taxonomy = tax.build_dag()
    k2_agg = Kraken2(
        taxonomy=taxonomy,
        report_folder=kraken_dir,
        output_file=output_csv,
        pattern=pattern,
    )
    k2_agg.aggregate_reports()
