from parsers.fastq_stats import FqLength
import click
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

'''
Parse length of the reads from a given fastq file and save to json format
'''
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--fq", "fastq", required=True, help="Fastq file (supports .gz files)"
)
@click.option(
    "--json", "json", required=True, help="File name for json formatted output file"
)
def read_length_parser(fastq: str, json: str) -> None:
    flp = FqLength(fq = fastq, out_file=json)
    flp.read_length()