import click
from statter.parsers.sourmash_kreport_parser import SourmashMetagenome

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--report_dir",
    "report_dir",
    required=True,
    help="Directory containing sourmash metagenome classification reports (see `sourmash tax metagenome -h`)",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
)
@click.option(
    "--pattern",
    "pattern",
    default="*.txt",
    help="Report naming pattern",
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
def sourmash_report_aggregator(report_dir: str, pattern: str, output_csv: str):
    """Aggregate `sourmash tax metagenome` classification reports into a single CSV"""
    parser = SourmashMetagenome(report_dir=report_dir,output_file=output_csv,pattern=pattern)
    parser.aggregate_reports()