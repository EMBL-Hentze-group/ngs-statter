import json
import logging
from pathlib import Path
from typing import List, Optional, Set

import pandas as pd

from statter.json_models.sample_stats import AllStats, StarStats

logger = logging.getLogger(__name__)


class SampleStats:
    def __init__(
        self,
        sample: str,
        raw: str,
        first_trim: str,
        align: str,
        out: str,
        second_trim: Optional[str] = None,
        rRNA_free: Optional[str] = None,
        rRNA_mapped: Optional[str] = None,
        dedup: Optional[str] = None,
        kraken2: Optional[str] = None,
    ) -> None:
        self.out = out
        self._all_stats = AllStats(
            sample_name=sample,
            raw_reads=self._seqkit_reader(raw),
            first_trim=self._seqkit_reader(first_trim),
        )  # type: ignore
        if second_trim:
            self._all_stats.second_trim = self._seqkit_reader(second_trim)
        if rRNA_free:
            self._all_stats.rRNA_free = self._seqkit_reader(rRNA_free)
        if rRNA_mapped:
            self._all_stats.rRNA_mapped = self._seqkit_reader(rRNA_mapped)
        self._align_stats: StarStats = self._json_reader(align)
        self._update_aln_stats()
        if dedup:
            self._dedup_stats: StarStats = self._json_reader(dedup)
            self._update_dedup_stats()
        self._classified: Optional[int] = None
        self._unclassified: Optional[int] = None
        if kraken2:
            self._parse_kraken2_report(kraken2)
            self._all_stats.classified = self._classified
            self._all_stats.unclassified = self._unclassified

    @staticmethod
    def _seqkit_reader(fname: str) -> int:
        """_seqkit_reader Helper function
        return num_seqs column. If the reads are paired end, return the max. value

        Args:
            fname: str, tsv formatted, from  seqkit stats -aT

        Returns:
            int, read count
        """
        return pd.read_csv(fname, sep="\t", usecols=["num_seqs"]).max().item()

    @staticmethod
    def _json_reader(fname: str) -> StarStats:
        """_json_reader Helper function
        Return alignment and deduplicated bam statistics

        Args:
            fname: Bam file stats in json format

        Returns:
            StarStats, alignment and deduplicated bam statistics
        """
        with open(fname, "r") as jh:
            jdata = json.load(jh)
        stats: StarStats = StarStats.model_validate(jdata)
        return stats

    def _update_aln_stats(self) -> None:
        """_update_aln_stats
        Update alignment statistics to the main AllStats object
        """
        self._all_stats = self._all_stats.model_copy(
            update=self._align_stats.model_dump()
        )

    def _update_dedup_stats(self) -> None:
        """_update_dedup_stats
        Update deduplicated alignment statistics to the main AllStats object
        """
        self._all_stats.umi_mapped_total = self._dedup_stats.mapped_total
        self._all_stats.umi_mapped_uniquely_mapped_reads = (
            self._dedup_stats.mapped_uniquely_mapped_reads
        )
        self._all_stats.umi_mapped_multimapped_reads = (
            self._dedup_stats.mapped_multimapped_reads
        )

    def _parse_kraken2_report(self, kraken2: str) -> None:
        """_parse_kraken2_report Helper function
        Parse kraken2 report to get classified and unclassified read counts
        Kranen2 report format: https://github.com/DerrickWood/kraken2/wiki/Manual
        """
        report: pd.DataFrame = pd.read_csv(
            kraken2,
            sep="\t",
            header=None,
            names=[
                "percent",
                "num_reads_clade",
                "num_reads_taxon",
                "rank",
                "tax_id",
                "name",
            ],
        )
        try:
            self._classified = int(
                report.loc[report["rank"] == "R", "num_reads_clade"].item()
            )
        except ValueError:
            logging.warning(
                f"Could not find classified reads in kraken2 report: {kraken2}!"
            )
            pass
        try:
            self._unclassified = int(
                report.loc[report["rank"] == "U", "num_reads_clade"].item()
            )
        except ValueError:
            logging.warning(
                f"Could not find unclassified reads in kraken2 report: {kraken2}!"
            )
            pass

    def to_json(self) -> None:
        """to_json
        Write the complete sample statistics to a json file
        """
        if Path(self.out).exists():
            logger.warning(f"Overwriting existing stats file: {self.out}")
        with open(self.out, "w") as _jh:
            _jh.write(self._all_stats.model_dump_json(indent=1))


def compile_all_sample_stats(output: str, files: List[str]) -> None:
    """compile_all_sample_stats
    Collect statistics from multiple samples and compile into a single tsv file

    Args:
        output: str, output tsv file name
        files: list of str, list of json files to compile
    """
    stats_list = []
    for f in sorted(files):
        stats = AllStats.model_validate_json(Path(f).read_text())
        stats_list.append(stats.model_dump(by_alias=True))
    samples: Set[str] = set([s["Sample name"] for s in stats_list])
    if len(samples) != len(stats_list):
        raise ValueError(
            f"Duplicate sample names found: {len(stats_list)} files but only {len(samples)} unique samples!"
        )
    stats_df: pd.DataFrame = (
        pd.DataFrame(stats_list).set_index("Sample name").sort_index().transpose()
    )
    if Path(output).exists():
        logger.warning(f"Overwriting existing stats file: {output}")
    # replace NaN with pd.NA for proper tsv output
    stats_df.to_csv(output, sep="\t", index=True, na_rep="<NA>")
