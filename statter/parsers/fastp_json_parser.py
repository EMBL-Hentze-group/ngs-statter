from ast import Set
import logging
import json
from collections import defaultdict
from pathlib import Path
import re
from typing import List, Dict, Set

"""
parse fastp json outputs and generate trimming stats per sample report file
"""

logger = logging.getLogger(__file__)


class Fastp:
    """
    Fastp reports parser
    """

    def __init__(
        self,
        trim_dir: str,
        out_csv: str,
        first_trim_pattern: str = "*.first_trim.json",
        second_trim_pattern: str = "*.second_trim.json",
    ) -> None:
        self.trim_dir = trim_dir
        self.out_csv = out_csv
        if first_trim_pattern == second_trim_pattern:
            raise RuntimeError(
                f"File naming patterns first trim pattern '{first_trim_pattern}' and second trim pattern '{second_trim_pattern}' cannot be the same!"
            )
        self.first_trim_pattern = first_trim_pattern
        self.second_trim_pattern = second_trim_pattern
        self._trim1_files = self._collect_sample_files(pattern=first_trim_pattern)
        self._trim2_files = self._collect_sample_files(pattern=second_trim_pattern)
        if len(self._trim1_files) != len(self._trim2_files):
            logger.warning(
                f"Unequal number of files with first and second trim data! # first trim files: {len(self._first_trim_files)}, # second trim files: {len(self._second_trim_files)}"
            )
        # dictionary for trimming stats
        self.trimming_stats = defaultdict(dict)
        # set of all used keys
        self.stats_keys = set()

    def _list_json_files(self, pattern: str) -> List[Path]:
        """
        helper function
        """
        files = list(Path(self.trim_dir).glob(pattern))
        if len(files) == 0:
            raise RuntimeError(
                f"Cannot parse files with pattern '{pattern}' from {self.trim_dir}"
            )
        logger.info(f"Found {len(files)} with pattern '{pattern}' in {self.trim_dir} ")
        return files

    def _collect_sample_files(self, pattern: str) -> Dict[str, str]:
        """
        helper function
        glob and list files using the given pattern and
        collect fastp trim stat json file per sample
        """
        sample_files = {}
        pattern_clean = re.sub(r"\*{1,}", "", pattern) + ".*$"
        for f in Path(self.trim_dir).glob(pattern=pattern):
            sample = re.sub(pattern_clean, "", f.name)
            sample_files[sample] = str(f)
        if len(sample_files) == 0:
            raise RuntimeError(
                f"Cannot parse files with pattern '{pattern}' from {self.trim_dir}"
            )
        logger.info(
            f"Found {len(sample_files)} with pattern '{pattern}' in {self.trim_dir} "
        )
        return sample_files

    def collect_trimming_stats(self):
        """
        collect trimming data from accumulated first and second trim files
        """
        common = set(self._trim1_files.keys()) & set(self._trim2_files.keys())
        # do not use these keys!
        omit_keys = set(["passed_filter_reads"])
        if len(common) == 0:
            raise RuntimeError(
                f"Cannot find any common samples between patterns '{self.first_trim_pattern}' and '{self.second_trim_pattern}' in folder {self.trim_dir}!"
            )
        diff1 = common - set(self._trim1_files.keys())
        if len(diff1) > 0:
            diff1_samples = ", ".join(diff1)
            logging.warning(
                f"Skipping samples {diff1_samples}, cannot find corresponding second trimming files with pattern '{self.second_trim_pattern}' in {self.trim_dir}"
            )
        diff2 = common - set(self._trim2_files.keys())
        if len(diff2) > 0:
            diff2_samples = ", ".join(diff2)
            logging.warning(
                f"Skipping samples {diff2_samples}, cannot find corresponding first trimming files with pattern '{self.first_trim_pattern}' in {self.trim_dir}"
            )
        for sample in sorted(common):
            first_trim = self._load_json(self._trim1_files[sample])
            second_trim = self._load_json(self._trim2_files[sample])
            # Raw data
            raw_counts, trimmed_counts = 0, 0
            try:
                raw_counts = first_trim["summary"]["before_filtering"]["total_reads"]
            except KeyError:
                logger.warning(
                    "Cannot find key ['summary']['before_filtering']['total_reads']! replacing value with 0!"
                )
                pass
            # Final trimmed data
            try:
                trimmed_counts = second_trim["filtering_result"]["passed_filter_reads"]
            except KeyError:
                logger.warning(
                    "Cannot find key ['filtering_result']['passed_filter_reads']! replacing value with 0!"
                )
                pass
            self.trimming_stats[sample]["raw_counts"] = raw_counts
            self.trimming_stats[sample]["trimmed_counts"] = trimmed_counts
            # add in between data
            self._add_data(first_trim["filtering_result"], sample, omit_keys)
            self._add_data(second_trim["filtering_result"], sample, omit_keys)
        self._write_stats()

    def _load_json(self, fname: str) -> Dict:
        """
        helper function
        load json
        """
        with open(fname, "r") as jh:
            trim_dat = json.load(jh)
        return trim_dat

    def _add_data(self, json_dict: Dict, sample: str, omit_keys: Set[str]) -> None:
        """
        Helper function add data to respective sample key
        """
        for k, v in json_dict.items():
            if k in omit_keys or v == 0:
                continue
            self.stats_keys.add(k)
            try:
                self.trimming_stats[sample][k] += v
            except KeyError:
                self.trimming_stats[sample][k] = v

    def _write_stats(self) -> None:
        """
        Helper function
        write trimming stats to file
        """
        if Path(self.out_csv).exists():
            logger.warning(f"Re-writing file {self.out_csv}")
        headers = ["raw_counts", "trimmed_counts"] + sorted(self.stats_keys)
        with open(self.out_csv, "w") as wh:
            header_str = "\t".join(headers)
            wh.write(f"samples\t{header_str}\n")
            for s in sorted(self.trimming_stats.keys()):
                sdat = ["0"] * len(headers)
                for i, h in enumerate(headers):
                    try:
                        sdat[i] = str(self.trimming_stats[s][h])
                    except KeyError:
                        continue
                sdat_str = "\t".join(sdat)
                wh.write(f"{s}\t{sdat_str}\n")
