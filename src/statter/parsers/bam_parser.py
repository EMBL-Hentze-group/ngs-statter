import json
import logging
from pathlib import Path
from typing import Any, Dict

from statter.statter import alignment_stats, gene_type_read_dist, star_bam_stats

from statter.json_models.sample_stats import BamStats, StarStats


class BamParser:

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
            "Reads for mapping": 0,
            "Mapped: Total": 0,
            "Mapped: Multimapped reads": 0,
            "Mapped: Uniquely mapped reads": 0,
            "Mapped: PCR duplicate reads": 0,
            "Mapped: Unique reads": 0,
            "Unmapped: Total": 0,
        }  # place holders for STAR alignment stats
        self._must_keys = set(
            [
                "Reads for mapping",
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
        - Unmapped: The number of reads that are not mapped to the reference genome.

        The alignment statistics are written to the specified output file in JSON format.
        """
        align_stats: dict[str, int] = alignment_stats(self.bam, self.min_q)
        self._to_json(BamStats(**align_stats).model_dump(), out_json)

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
