import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List
import logging

import pysam
from statter.parsers.gff_parser import Gene


class BamParser:
    """
    Parse bam file for data
    """

    def __init__(
        self, bam: str, min_q: int = 0, ignore_duplicate: bool = False
    ) -> None:
        """
        Arguments:
         bam: bam file path
         min_q: minimum alignment quality (int)
         ignore_duplicate: Flag to ignore PCR duplictes
        """
        self.bam = bam
        self.min_q = min_q
        self.ignore_duplicate = ignore_duplicate

    @staticmethod
    def _to_json(stat_data, out_file) -> None:
        """
        Helper function
        dump the given data in json format
        """
        if Path(out_file).exists():
            logging.warning(f"Over-writing file {out_file}")
        with open(out_file, "w") as jh:
            json.dump(stat_data, jh)

    def read_length_stats_per_gene_type(
        self, genes: Dict[str, List[Gene]], out_json: str
    ) -> None:
        """
        compute read length stats per gene
        """
        gene_read_lens = defaultdict(dict)
        with pysam.AlignmentFile(
            self.bam, mode="rb", require_index=True, duplicate_filehandle=True
        ) as bh:
            ref_chroms = set(bh.references)
            common_chroms = ref_chroms & set(genes.keys())
            if len(common_chroms) == 0:
                raise RuntimeError(
                    f"Cannot find any common chromosomes between {self.bam} and given gene data!"
                )
            for chrom in common_chroms:
                for gene in genes[chrom]:
                    for aln in bh.fetch(contig=chrom, start=gene.start, stop=gene.stop):
                        if (
                            aln.is_qcfail
                            or aln.is_secondary
                            or aln.is_supplementary
                            or aln.is_unmapped
                            or aln.mapping_quality < self.min_q
                        ) or (self.ignore_duplicate and aln.is_duplicate):
                            continue
                        strand = "-" if aln.is_reverse else "+"
                        if strand != gene.strand:
                            continue
                        # status = "dup" if aln.is_duplicate else "non_dup"
                        # key = (
                        #     gene.gene_id,
                        #     gene.name,
                        #     gene.gene_type,
                        #     status,
                        # )
                        try:
                            gene_read_lens[gene.gene_type][aln.query_length] += 1
                        except KeyError:
                            gene_read_lens[gene.gene_type][aln.query_length] = 1
        if len(gene_read_lens) == 0:
            raise RuntimeError(
                f"Cannot parse data form given bam file: {self.bam}! Check your input GFF3 and Bam file"
            )
        self._to_json(gene_read_lens, out_json)
