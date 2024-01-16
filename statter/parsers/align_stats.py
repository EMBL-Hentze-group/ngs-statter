import logging

"""
Compute alignment stats from raw fastq file to final genome aligned and genome unaligned reads
"""


class AlignmentStats:
    """
    compute alignment stats
    """

    def __init__(
        self,
        raw_fastq: str,
        trim_fastq: str,
        rrna_bam: str,
        rrna_unaligned_fastq: str,
        genome_bam: str,
        genome_unaligned_fastq: str,
    ) -> None:
        """
        raw_fastq: raw fastq file
        trim_fastq: trimmed fastq file
        rrna_bam: read alignment to rRNA/repbase
        rrna_unaligned_bam: reads not aligned to rRNA/repbase
        genome_bam: read alignment to genome
        genome_unaligned_fastq: reads not aligned to genome
        """
