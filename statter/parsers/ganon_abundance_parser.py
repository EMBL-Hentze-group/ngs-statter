import logging
from math import prod
from pathlib import Path
from collections import OrderedDict
from itertools import chain, product
import re

"""
ganon .tre file parser and result aggregator
"""
logger = logging.getLogger(__name__)


class GanonTreParser:
    """
    Parse and aggregate .tre abundance estimation results from ganon
    """

    def __init__(
        self, result_folder: str, output_file: str, pattern: str = "*.tre"
    ) -> None:
        self.result_folder = result_folder
        self.output_file = output_file
        self.pattern = pattern
        self.abundace = OrderedDict()
        self.sample_names = list()

    def aggregate_results(self) -> None:
        """
        Helper function
        load abundace estimation file
        Columns:
        rank : used column
        target -> taxonomic id
        lineage
        name : used column
        # unique : number of reads that matched exclusively to this target
        # shared
        # children
        # cumulative : used column
        % cumulative : used column
        """
        abundace_files = sorted(Path(self.result_folder).glob(self.pattern))
        if len(abundace_files) == 0:
            raise RuntimeError(
                f"Cannot find ganon report files with pattern {self.pattern} in {self.result_folder}"
            )
        for i, tre in enumerate(abundace_files):
            self.sample_names.append(re.sub(r"\..*$", "", tre.stem))
            j = 1
            with open(tre, "r") as _fh:
                for f in _fh:
                    fdat = f.strip().split("\t")
                    if len(fdat) < 9:
                        continue
                    uniq_id = (fdat[0], fdat[3], fdat[1])
                    try:
                        self.abundace[uniq_id][i] = (fdat[7], fdat[8])
                    except KeyError:
                        self.abundace[uniq_id] = [("0", "0")] * len(abundace_files)
                        self.abundace[uniq_id][i] = (fdat[7], fdat[8])
        self._write_results()

    def _write_results(self):
        """
        helper function
        write results
        """
        if Path(self.output_file).exists():
            logger.warning(f"Re-writing file {self.output_file}")
        sample_names = ["_".join(s) for s in product(self.sample_names, ["count", "%"])]
        header = "\t".join(["rank", "name", "taxonomy_id"] + sample_names)
        with open(self.output_file, "w") as wh:
            wh.write(header + "\n")
            for species, counts in self.abundace.items():
                row = "\t".join(list(species) + list(chain(*counts)))
                wh.write(row + "\n")
