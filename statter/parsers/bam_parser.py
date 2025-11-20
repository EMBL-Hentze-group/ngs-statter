import json
import logging
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Any, Dict, List

import pysam

from statter.json_models.sample_stats import StarStats
from statter.parsers.gff_parser import Gene
from statter.statter import alignment_stats, gene_type_read_dist, star_bam_stats


class BamParser:
    """
    Parse bam file for data.

    Attributes:
        bam (str): The bam file name. Must be co-ordinate sorted and indexed.
        min_q (int): Minimum alignment quality. Defaults to 0.
        ignore_duplicate (bool): Flag to ignore PCR duplicates (if given by `samtools markdup`). Defaults to False.

    Methods:
        __init__(self, bam: str, min_q: int = 0, ignore_duplicate: bool = False) -> None:
            Initialize the BamParser object.

        _to_json(stat_data: Any, out_file: str) -> None:
            Helper function to dump the given data in json format.

        read_length_stats_per_gene_type(self, genes: Dict[str, List[Gene]], out_json: str) -> None:
            Compute read length stats per gene.

        STAR_alignment_stats(self, out_json: str) -> None:
            Write STAR alignment stats to json.

        STAR_alignment_stats_rs(self, out_json: str) -> None:
            Write STAR alignment stats to json using the `star_bam_stats` function.

        alignment_stats(self, out_json: str) -> None:
            Write alignment stats to json using the `alignment_stats` function.

        _unmapped_type(self, ut_type: str) -> str:
            Helper function to return STAR unmapped type.

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
        self._map_stats: Dict[str, int] = {
            "Input reads": 0,
            "Mapped: Total": 0,
            "Mapped: Multimapped reads": 0,
            "Mapped: Uniquely mapped reads": 0,
            "Mapped: PCR duplicate reads": 0,
            "Mapped: Unique reads": 0,
            "Unmapped: Total": 0,
        }  # place holders for STAR alignment stats
        self._must_keys = set(
            [
                "Input reads",
                "Mapped: Total",
                "Mapped: Multimapped reads",
                "Mapped: Uniquely mapped reads",
            ]
        )  # these are mandatory keys

    @staticmethod
    def _to_json(stat_data: Any, out_file: str) -> None:
        """
        Helper function
        dump the given data in json format
        """
        if Path(out_file).exists():
            logging.warning(f"Over-writing file {out_file}")
        with open(out_file, "w") as jh:
            json.dump(stat_data, jh, indent=1)

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

    def read_length_stats_per_gene_type_rs(
        self,
        gff3_file: str,
        out_json: str,
        features: list = ["gene", "tRNA"],
        gene_id: str = "gene_id",
        gene_name: str = "gene_name",
        gene_type: str = "gene_type",
    ) -> None:
        gene_type_dist: dict[str, dict[int, int]] = gene_type_read_dist(
            self.bam,
            self.min_q,
            self.ignore_duplicate,
            gff3_file,
            features,
            gene_id,
            gene_name,
            gene_type,
        )
        if len(gene_type_dist) == 0:
            raise RuntimeError(
                f"Cannot find reads in {self.bam} aligned to features {features} in {gff3_file}"
            )
        self._to_json(gene_type_dist, out_json)

    def STAR_alignment_stats(self, out_json: str) -> None:
        """STAR_alignment_stats
        Parse bam file from STAR aligner for alignment stats

        Args:
            out_json: json file name to write the alignment stats
        """
        with pysam.AlignmentFile(self.bam, "rb") as bh:
            for aln in bh:
                if aln.is_paired and not aln.is_read1:
                    continue
                if aln.is_unmapped:
                    self._map_stats["Unmapped: Total"] += 1
                    self._map_stats["Input reads"] += 1
                    try:
                        ut_tag = aln.get_tag("uT")
                    except KeyError:
                        logging.warning(
                            "Cannot find 'uT' tag in the bam file! Check if the bam file was generated using STAR aligner."
                        )
                        continue
                    unmapped_tag = self._unmapped_type(str(ut_tag))
                    try:
                        self._map_stats[unmapped_tag] += 1
                    except KeyError:
                        self._map_stats[unmapped_tag] = 1
                    continue
                if (
                    aln.is_qcfail
                    or aln.is_secondary
                    or aln.is_supplementary
                    or aln.mapping_quality < self.min_q
                ):
                    continue
                self._map_stats["Input reads"] += 1
                self._map_stats["Mapped: Total"] += 1
                if aln.is_duplicate:
                    self._map_stats["Mapped: PCR duplicate reads"] += 1
                else:
                    self._map_stats["Mapped: Unique reads"] += 1
                try:
                    n_aln = aln.get_tag("NH")
                except KeyError:
                    n_aln = 0
                    pass
                if n_aln == 1:
                    self._map_stats["Mapped: Uniquely mapped reads"] += 1
                elif n_aln > 1:  # type: ignore
                    self._map_stats["Mapped: Multimapped reads"] += 1
        self._to_json(StarStats(**self._map_stats).model_dump(), out_json)

    def STAR_alignment_stats_rs(self, out_json: str) -> None:
        """STAR_alignment_stats_rs
        Parse bam file from STAR aligner for alignment stats using rust function
        Args:
            out_json: json file name to write the alignment stats

        Raises:
            RuntimeError: If alignment statistics cannot be parsed from the BAM file.
            RuntimeError: If any required values are missing from the alignment statistics.
        """
        self._map_stats = star_bam_stats(self.bam, self.min_q)
        if len(self._map_stats) == 0:
            raise RuntimeError(
                f"Cannot parse alignment stats from {self.bam}! Check your input file"
            )
        diff_keys = self._must_keys - set(self._map_stats.keys())
        if len(diff_keys) > 0:
            missing_keys = ", ".join(diff_keys)
            raise RuntimeError(
                f"Cannot find the following values: {missing_keys} from {self.bam}. Check your input bam file!"
            )
        self._to_json(StarStats(**self._map_stats).model_dump(), out_json)

    def alignment_stats_rs(self, out_json: str) -> None:
        """
        Compute alignment statistics for a BAM file.

        Args:
            out_json (str): The output file name for the alignment statistics in JSON format.

        Raises:
            RuntimeError: If alignment statistics cannot be parsed from the BAM file.
            RuntimeError: If any required values are missing from the alignment statistics.
            RuntimeError: If the specified values cannot be found in the BAM file.

        Returns:
            None

        The alignment statistics include the following metrics:
        - Input reads: The total number of reads in the BAM file.
        - Mapped: The number of reads that are mapped to the reference genome.
        - Mapped %: The percentage of reads that are mapped to the reference genome.
        - Unmapped: The number of reads that are not mapped to the reference genome.
        - Unmapped %: The percentage of reads that are not mapped to the reference genome.

        The alignment statistics are written to the specified output file in JSON format.
        """
        align_stats: dict[str, int] = alignment_stats(self.bam, self.min_q)
        if len(align_stats) == 0:
            raise RuntimeError(
                f"Cannot parse alignment stats from {self.bam}! Check your input file"
            )
        map_keys = set(["Input reads", "Mapped", "Unmapped"])
        diff_keys = map_keys - set(align_stats.keys())
        if len(diff_keys) > 0:
            missing_keys = ", ".join(diff_keys)
            raise RuntimeError(
                f"Cannot find the following values: {missing_keys} from {self.bam}. Check your input bam file!"
            )
        out_stats: OrderedDict[str, int | float] = OrderedDict()
        out_stats["Input reads"] = align_stats["Input reads"]
        # Mapped
        out_stats["Mapped"] = align_stats["Mapped"]
        out_stats["Mapped %"] = round(
            align_stats["Mapped"] * 100 / align_stats["Input reads"], 3
        )
        # Unmapped
        out_stats["Unmapped"] = align_stats["Unmapped"]
        out_stats["Unmapped %"] = round(
            align_stats["Unmapped"] * 100 / align_stats["Input reads"], 3
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
            "3": "Unmapped: mapped to too many loci",
            "4": "Unmapped: paired-end mate",
        }
        try:
            unmapped = unmapped_type[ut_type]
        except KeyError:
            unmapped = "Unknown"
        return unmapped
