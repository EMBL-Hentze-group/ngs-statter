from typing import List

import rich_click as click

from statter.parsers.bam_parser import BamParser
from statter.plotters.gene_type_read_length_distribution import GeneTypePlot

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def genetype() -> None:
    """
    Utilities for processing gene type read length distributions
    """


@genetype.command(
    "parse-gene-type-read-length",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
)
@click.option(
    "--gff3",
    "gff3",
    required=True,
    help="Gene annotation in GFF3 format (supports .gz files)",
)
@click.option(
    "--gene-features",
    "gene_features",
    default=["gene", "tRNA"],
    help="Features to parse from GFF3 file, MUST be supplied as: '--gene-features gene --gene-features tRNA'...",
    multiple=True,
    show_default=True,
    type=str,
)
@click.option(
    "--gene-id",
    "gene_id",
    default="gene_id",
    help="Gene id attribute in GFF3 attribute column",
    show_default=True,
    type=str,
)
@click.option(
    "--gene-name",
    "gene_name",
    default="gene_name",
    help="Gene name attribute in GFF3 attribute column",
    show_default=True,
    type=str,
)
@click.option(
    "--gene-type",
    "gene_type",
    default="gene_type",
    help="Gene type attribute in GFF3 attribute column",
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
    "--min-q",
    "min_q",
    default=0,
    help="Minimum alignment quality",
    show_default=True,
    type=int,
)
@click.option(
    "--ignore-duplicate",
    "ignore_duplicate",
    is_flag=True,
    show_default=True,
    default=False,
    help="Flag to ignore PCR duplicate reads (after samtools markdup)",
)
@click.option(
    "--out-json",
    "out_json",
    required=True,
    help="Output file to write json formatted data",
)
def gene_type_read_length_stats(
    gff3: str,
    bam: str,
    out_json: str,
    gene_features: List[str],
    gene_id: str,
    gene_name: str,
    gene_type: str,
    min_q: int = 0,
    ignore_duplicate: bool = False,
) -> None:
    """
    Compute gene type read length distributions\b

    Given a gff3 file, a list of gene types to parse, and a bam file
    calculate the length distribution of reads aligning to these gene types
    and write this info in json format
    """
    bam_parser = BamParser(bam=bam, min_q=min_q, ignore_duplicate=ignore_duplicate)
    bam_parser.read_length_stats_per_gene_type_rs(
        gff3_file=gff3,
        out_json=out_json,
        features=gene_features,
        gene_id=gene_id,
        gene_name=gene_name,
        gene_type=gene_type,
    )


@genetype.command(
    "plot-gene-type-read-length",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
)
@click.argument(
    "json",
    nargs=-1,
    required=False,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    # help ="List of json formatted read length files (see parse_gene_type_read_length -h)",
)
@click.option(
    "--json-dir",
    "json_dir",
    required=False,
    help="Directory containing json formatted read length files (see statter parse-gene-type-read-length  -h)",
)
@click.option(
    "--pattern",
    "pattern",
    default="*.json",
    help="File naming pattern",
    show_default=True,
    type=str,
)
@click.option(
    "--nsamples",
    "nsamples",
    default=1,
    help="Number of replicates per group",
    show_default=True,
    type=int,
)
@click.option("--html", "output_html", required=True, help="Plot output file name")
def gene_type_length_plotter(
    json_dir: str, json: list[str], output_html: str, pattern: str, nsamples=1
) -> None:
    """
    Plot gene type read length distributions\b

    Read in json formatted gene type read length stats and output html plot.
    Either a directory containing json files (using --json-dir option) or a list of json files as space separated arguments MUST be provided.
    """
    if json_dir is None and len(json) == 0:
        raise click.UsageError(
            "Either --json-dir option or json formatted files as space separated arguments must be provided"
        )
    if json_dir is not None and len(json) > 0:
        raise click.UsageError(
            "Either --json-dir option or json formatted files as space separated arguments must be provided, not both"
        )
    gtl = GeneTypePlot(
        json_folder=json_dir,
        json_files=json,
        output_file=output_html,
        pattern=pattern,
        nsamples=nsamples,
    )
    gtl.plot()
