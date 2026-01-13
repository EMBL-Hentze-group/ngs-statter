from typing import List

import click

from statter.parsers.bam_parser import BamParser
from statter.plotters.gene_type_read_length_distribution import GeneTypePlot

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--bam",
    "bam",
    required=True,
    help="Alignments BAM file (MUST be co-ordinate sorted and indexed)",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
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
    type=click.Path(exists=False),
)
def alignment_stats_STAR(bam: str, min_q: int, out_json: str) -> None:
    bam_parser = BamParser(bam=bam, min_q=min_q, ignore_duplicate=False)
    bam_parser.STAR_alignment_stats_rs(out_json=out_json)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--bam",
    "bam",
    required=True,
    help="Alignments BAM file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
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
    type=click.Path(exists=False),
)
def basic_alignment_stats(bam: str, min_q: int, out_json: str) -> None:
    """
    Compute basic alignment statistics for a BAM file.

    \b
    Raises:
        RuntimeError: If alignment statistics cannot be parsed from the BAM file.
        RuntimeError: If any required values are missing from the alignment statistics.
        RuntimeError: If the specified values cannot be found in the BAM file.

    \b
    The basic alignment statistics include the following metrics:
    - Input reads: The total number of reads in the BAM file.
    - Mapped: The number of reads that are mapped to the reference genome.
    - Mapped %: The percentage of reads that are mapped to the reference genome.
    - Unmapped: The number of reads that are not mapped to the reference genome.
    - Unmapped %: The percentage of reads that are not mapped to the reference genome.

    The alignment statistics are written to the specified output JSON file.
    """
    bam_parser = BamParser(bam=bam, min_q=min_q, ignore_duplicate=False)
    bam_parser.alignment_stats_rs(out_json=out_json)


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
    help="Features to parse from GFF3 file, MUST be supplied as: '--gene_features gene --gene_features tRNA'...",
    multiple=True,
    show_default=True,
    type=str,
)
@click.option(
    "--gene_id",
    "gene_id",
    default="gene_id",
    help="Gene id attribute in GFF3 attribute column",
    show_default=True,
    type=str,
)
@click.option(
    "--gene_name",
    "gene_name",
    default="gene_name",
    help="Gene name attribute in GFF3 attribute column",
    show_default=True,
    type=str,
)
@click.option(
    "--gene_type",
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
    "--min_q",
    "min_q",
    default=0,
    help="Minimum alignment quality",
    show_default=True,
    type=int,
)
@click.option(
    "--ignore_duplicate",
    "ignore_duplicate",
    is_flag=True,
    show_default=True,
    default=False,
    help="Flag to ignore PCR duplicate reads (after samtools markdup)",
)
@click.option(
    "--out_json",
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
    Given a gff3 file, a list of features to parse, and a bam file
    calculate the length distribution of reads aligning to these features
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


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--json_dir",
    "json_dir",
    required=True,
    help="Directory containing json formatted read_length files (see gene_type_read_length_parser -h)",
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
    json_dir: str, output_html: str, pattern: str, nsamples=1
) -> None:
    """
    Read in json formatted gene type read length stats and output html plot
    """
    gtl = GeneTypePlot(
        json_folder=json_dir,
        output_file=output_html,
        pattern=pattern,
        nsamples=nsamples,
    )
    gtl.plot()
