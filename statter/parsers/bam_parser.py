import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List
import logging

import pysam
from gff_parser import Gene


class BamParser:
    """
    Parse bam file for data
    """

    def __init__(self, bam: str, min_q: int = 0) -> None:
        self.bam = bam
        self.min_q = min_q

    @staticmethod
    def _to_json(stat_data, out_file) -> None:
        """
        Helper function
        dump the given data in json format
        """
        if Path(out_file).exists():
            logging.warning(f"Over-writing file {out_file}")
        with open(out_file, "w") as jh:
            json.dump(stat_data, out_file)

    def read_length_stats_per_gene(
        self, genes: Dict[str, List[Gene]], out_json: str
    ) -> None:
        """
        compute read length stats per gene
        """
        gene_read_lens = defaultdict(Dict)
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
                        ):
                            continue
                        strand = "-" if aln.is_reverse else "+"
                        if strand != gene.strand:
                            continue
                        status = "dup" if aln.is_duplicate else "non_dup"
                        key = (
                            gene.gene_id,
                            gene.name,
                            gene.gene_type,
                            status,
                        )
                        try:
                            gene_read_lens[key][aln.query_length] += 1
                        except KeyError:
                            gene_read_lens[key][aln.query_length] = 1
        if len(gene_read_lens) == 0:
            raise RuntimeError(
                f"Cannot parse data form given bam file: {self.bam}! Check your input GFF3 and Bam file"
            )
        self._to_json(gene_read_lens, out_json)
