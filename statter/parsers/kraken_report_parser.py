import logging
import re
from itertools import chain, product
from pathlib import Path
from typing import List
from networkx import DiGraph, dfs_postorder_nodes

"""
aggregate per sample Kraken2 reports and generate one single file
"""

logger = logging.getLogger(__file__)


class Kraken2:
    def __init__(
        self,
        taxonomy: DiGraph,
        report_folder: str,
        output_file: str,
        pattern: str = "_report.txt",
    ) -> None:
        self.taxonomy = taxonomy
        self.output_file = output_file
        self.pattern = pattern
        self._kraken_files = sorted(Path(report_folder).glob(pattern))
        if len(self._kraken_files) == 0:
            raise RuntimeError(
                f"Cannot find Kraken2 report files with pattern {pattern} in folder {report_folder}"
            )
        logger.info(
            f"Found {len(self._kraken_files)} kraken2 report files in {report_folder}"
        )
        self._kraken_report = {}
        self._sample_names = []

    def aggregate_reports(self):
        """
        Aggregate Kraken2 reports into one dictionary
        Kraken report fields:
            1. Percentage of fragments covered by the clade rooted at this taxon
            2. Number of fragments covered by the clade rooted at this taxon
            3. Number of fragments assigned directly to this taxon
           *4. Number of minimizers in read data associated with this taxon (new)
           *5. An estimate of the number of distinct minimizers in read data associated with this taxon (new)
            6. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
            7. NCBI taxonomic ID number
            8. Indented scientific name
            * = Optional
        source: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
        """
        for i, kr in enumerate(self._kraken_files):
            self._sample_names.append(re.sub(r"\..*$", "", kr.stem))
            with open(kr, "r") as _kh:
                for l in _kh:
                    ldat = l.strip().split("\t")
                    try:
                        self._kraken_report[ldat[-2]][i] = (ldat[1], ldat[0], ldat[2])
                    except KeyError:
                        self._kraken_report[ldat[-2]] = [("0", "0", "0")] * len(
                            self._kraken_files
                        )
                        self._kraken_report[ldat[-2]][i] = (ldat[1], ldat[0], ldat[2])
        if len(self._kraken_report) == 0:
            kfiles = ", ".join([str(f) for f in self._kraken_files])
            raise RuntimeError(
                f"Cannot parse organism specific read classification report from file :{kfiles}. Check the input format!"
            )
        logger.info(
            f"Parsed {len(self._kraken_report)} taxonomy ranks from {len(self._kraken_files)} files"
        )
        self._write_report()

    def _write_report(self):
        """
        Helper function
        Write aggregated report to file
        """
        if Path(self.output_file).exists():
            logger.warning(f"Re-writing file {self.output_file}")
        snames = [
            "_".join(s)
            for s in product(
                self._sample_names, ["cumulative", "cumulative-%", "single"]
            )
        ]
        header = "\t".join(["order_id", "taxonomy_id", "name", "rank"] + snames)
        order_count = 0  # column to sort back into taxonomy order
        with open(self.output_file, "w") as wh:
            wh.write(f"{header}\n")
            try:
                unclassified = self._kraken_report["0"]
                out_str = "\t".join(
                    [f"{order_count}", "0", "unclassified", "no rank"]
                    + list(chain(*unclassified))
                )
                wh.write(f"{out_str}\n")
            except KeyError:
                logger.warning(
                    f"Cannot find 'unclassified' read data in Kraken reports!"
                )
                pass
            order_count += 1
            for idx in dfs_postorder_nodes(self.taxonomy):
                if idx not in self._kraken_report:
                    continue
                out_str = "\t".join(
                    [
                        f"{order_count}",
                        idx,
                        self.taxonomy.nodes[idx]["name"],
                        self.taxonomy.nodes[idx]["rank"],
                    ]
                    + list(chain(*self._kraken_report[idx]))
                )
                wh.write(f"{out_str}\n")
                order_count += 1
