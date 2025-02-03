import logging
import re
from itertools import chain, product
from pathlib import Path
from typing import List

""" 
aggregate per sample output from `sourmash tax metagenome` and collate the reports into a single file
"""

logger = logging.getLogger(__file__)

class SourmashMetagenome:
    def __init__(self, report_dir:str, output_file: str, pattern:str = "_report.txt") -> None:
        # self.report_dir = report_dir
        self.output_file = output_file
        # self.pattern = pattern
        self._reports = sorted(Path(report_dir).glob(pattern))
        if len(self._reports)==0:
            raise RuntimeError(
                f"Cannot find sourmash metagenome report files with pattern {pattern} in folder {report_dir}"
            )
        logger.info(
            f"Found {len(self._reports)} sourmash metagenome  report files in {report_dir}"
        )
        self._metagenome_report = {}
        self._sample_names: List[str] = []
    
    def aggregate_reports(self):
        """
        Aggregate sourmash metagenome kraken style reports into one dictionary
        Report fields:
            1. Percentage of fragments covered by the clade rooted at this taxon
            2. Number of fragments covered by the clade rooted at this taxon
            3. Number of fragments assigned directly to this taxon
            4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank. E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
            5. *NCBI taxonomic ID number
            6. Indented scientific name
            * = Optional
        source: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
        """
        for i,report in enumerate(self._reports):
            self._sample_names.append(re.sub(r"\..*$", "", report.stem))
            with open(report,"r") as _rh:
                for l in _rh:
                    ldat = l.strip().split("\t")
                    try:
                        self._metagenome_report[(ldat[5],ldat[3])][i] = (ldat[1], ldat[0], ldat[2])
                    except KeyError:
                        self._metagenome_report[(ldat[5],ldat[3])] = [("0", "0", "0")] * len(
                            self._reports
                        )
                        self._metagenome_report[(ldat[5],ldat[3])][i] = (ldat[1], ldat[0], ldat[2])
        if len(self._metagenome_report) == 0:
            kfiles = ", ".join([str(f) for f in self._reports])
            raise RuntimeError(
                f"Cannot parse organism specific read classification report from file(s) :{kfiles}. Check the input format!"
            )
        logger.info(
            f"Parsed {len(self._metagenome_report)} taxonomy ranks from {len(self._reports)} files"
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
                self._sample_names, ["cumulative", "cumulative-%", "individual"]
            )
        ]
        header = "\t".join(["name","rank"]+snames)
        with open(self.output_file, "w") as wh:
            wh.write(header + "\n")
            for tax, dat in self._metagenome_report.items():
                out = "\t".join(list(tax)+list(chain(*dat)))
                wh.write(out+"\n")