import click
from typing import List
from parsers.bam_parser import BamParser
from parsers.gff_parser import GFF3

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--gff3",
    "gff3",
    required=True,
    help="Gene annotation in GFF3 format (supports .gz files)",
)
@click.option(
    "--gene_features",
    "gene_features",
    default=["gene", "tRNA"],
    help="Features to parse from GFF3 file",
    show_default=True,
    type=str,
)
@click.option(
    "--bam",
    "bam",
    required=True,
    help="Alignments BAM file (MUST be co-ordinate sorted and indexed)",
)
@click.option(
    "--min_q",
    "min_q",
    default=0,
    help="Minimum alignment quality",
    show_default=True,
    type=int,
)
@click.option(
    "--out_json",
    "out_json",
    required=True,
    help="Output file to write json formatted data",
)
def gene_read_len_stats(
    gff3: str,
    bam: str,
    out_json: str,
    gene_features: List[str] = ["gene", "tRNA"],
    min_q: int = 0,
) -> None:
    """
    Given a gff3 file, a list of features to parse, and a bam file
    calculate the length distribution of reads aligning to these features
    and write this info in json format
    """
    gff3_parse = GFF3(gff3=gff3, features=gene_features)
    gff3_genes = gff3_parse.parse_gene_info()
    bam_parser = BamParser(bam=bam, min_q=min_q)
    bam_parser.read_length_stats_per_gene(genes=gff3_genes, out_json=out_json)
