import json
from collections import defaultdict, OrderedDict
from pathlib import Path
from typing import Dict, List
import logging

import pysam
from statter.parsers.gff_parser import Gene
from statter.statter import star_bam_stats


class BamParser:
    """
    Parse bam file for data
    """

    def __init__(
        self, bam: str, min_q: int = 0, ignore_duplicate: bool = False
    ) -> None:
        """__init__ BamPaser init

        Args:
            bam: bam file name. MUST BE co-ordinate sorted and indexed
            min_q: Minimum alignment qualty. Defaults to 0.
            ignore_duplicate: Flag to ignore PCR duplicates (if given by ``samtools markdup``). Defaults to False.
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

    def STAR_alignment_stats(self, out_json: str) -> None:
        """
        Write STAR alignment stats to json
        """
        map_stats = {"Unique": 0, "PCR duplicate": 0, "Total aligned": 0}
        with pysam.AlignmentFile(self.bam, "rb") as bh:
            for aln in bh:
                if aln.is_unmapped:
                    unmapped_tag = self._unmapped_type(str(aln.get_tag("uT")))
                    try:
                        map_stats[unmapped_tag] += 1
                    except KeyError:
                        map_stats[unmapped_tag] = 1
                else:
                    map_stats["Total aligned"] += 1
                    if aln.is_duplicate:
                        map_stats["PCR duplicate"] += 1
                    else:
                        map_stats["Unique"] += 1
        self._to_json(map_stats, out_json)

    def STAR_alignment_stats_rs(self, out_json: str) -> None:
        map_stats = star_bam_stats(self.bam, self.min_q)
        if len(map_stats) == 0:
            raise RuntimeError(
                f"Cannot parse alignment stats from {self.bam}! Check your input file"
            )
        map_keys = set(
            [
                "Unmapped: too short",
                "Unmapped: Total",
                "Mapped: Uniquely mapped reads",
                "Mapped: PCR duplicate reads",
                "Mapped: Unique reads",
                "Multimapping: mapped to too many loci",
                "Unmapped: no seed/windows",
                "Mapped: Total",
                "Input reads",
                "Mapped: Multimapped reads",
            ]
        )
        diff_keys = map_keys - set(map_stats.keys())
        if len(diff_keys) > 0:
            missing_keys = ", ".join(diff_keys)
            raise RuntimeError(
                f"Cannot find the following values: {missing_keys} from {self.bam}. Check your input bam file!"
            )
        out_stats: OrderedDict[str, int | float] = OrderedDict()
        out_stats["Input reads"] = map_stats["Input reads"]
        out_stats["Mapped: Total"] = map_stats["Mapped: Total"]
        # mapped
        out_stats["Mapped: Total %"] = round(
            float(map_stats["Mapped: Total"]) * 100 / map_stats["Input reads"], 3
        )
        # unique reads
        out_stats["Mapped: Unique reads"] = map_stats["Mapped: Unique reads"]
        out_stats["Mapped: Unique reads %"] = round(
            float(map_stats["Mapped: Unique reads"]) * 100 / map_stats["Input reads"], 3
        )
        # pcr duplicates
        out_stats["Mapped: PCR duplicate reads"] = map_stats[
            "Mapped: PCR duplicate reads"
        ]
        out_stats["Mapped: PCR duplicate reads %"] = round(
            float(map_stats["Mapped: PCR duplicate reads"])
            * 100
            / map_stats["Input reads"],
            3,
        )
        # uniquely mapped reads
        out_stats["Mapped: Uniquely mapped reads"] = map_stats[
            "Mapped: Uniquely mapped reads"
        ]
        out_stats["Mapped: Uniquely mapped reads %"] = round(
            float(map_stats["Mapped: Uniquely mapped reads"])
            * 100
            / map_stats["Input reads"],
            3,
        )
        # multimapped reads
        out_stats["Mapped: Multimapped reads"] = map_stats["Mapped: Multimapped reads"]
        out_stats["Mapped: Multimapped reads %"] = round(
            float(map_stats["Mapped: Multimapped reads"])
            * 100
            / map_stats["Input reads"],
            3,
        )
        # unmapped
        out_stats["Unmapped: Total"] = map_stats["Unmapped: Total"]
        out_stats["Unmapped: Total %"] = round(
            float(map_stats["Unmapped: Total"]) * 100 / map_stats["Input reads"], 3
        )
        # Multimapping: mapped to too many loci
        out_stats["Unmapped: mapped to too many loci"] = map_stats[
            "Multimapping: mapped to too many loci"
        ]
        out_stats["Unmapped: mapped to too many loci %"] = round(
            float(map_stats["Multimapping: mapped to too many loci"])
            * 100
            / map_stats["Input reads"],
            3,
        )
        # Unmapped: no seed/windows
        out_stats["Unmapped: no seed/windows"] = map_stats["Unmapped: no seed/windows"]
        out_stats["Unmapped: no seed/windows %"] = round(
            float(map_stats["Unmapped: no seed/windows"])
            * 100
            / map_stats["Input reads"],
            3,
        )
        # Unmapped: too short
        out_stats["Unmapped: too short"] = map_stats["Unmapped: too short"]
        out_stats["Unmapped: too short %"] = round(
            float(map_stats["Unmapped: too short"]) * 100 / map_stats["Input reads"],
            3,
        )
        self._to_json(out_stats, out_json)

    def _unmapped_type(self, ut_type: str) -> str:
        """
        Helper function
        Return STAR unmapped type:
        see documentation: https://raw.githubusercontent.com/alexdobin/STAR/master/doc/STARmanual.pdf
        """
        unmapped_type = {
            "0": "Unmapped: no seed/windows",
            "1": "Unmapped: too short",
            "2": "Unmapped: too many mismatches",
            "3": "Multimapping: mapped to too many loci",
            "4": "Unmapped: paired-end mate",
        }
        try:
            unmapped = unmapped_type[ut_type]
        except KeyError:
            unmapped = "Unknown"
        return unmapped
