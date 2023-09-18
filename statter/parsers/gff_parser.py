import gzip
import logging
import re
from dataclasses import dataclass
from functools import partial
from typing import Callable, Dict, List, Set

"""
    GFF3 parser
"""


@dataclass
class Gene:
    gene_id: str
    name: str
    gene_type: str
    chr: str
    start: int
    stop: int
    strand: str


class GFF3:
    """
    Parse GFF3 file and return gene info per chromosome
    """

    def __init__(self, gff3: str, features=["gene", "tRNA"]) -> None:
        self.gff3: str = gff3
        self._reader: Callable = self._parser()
        self._gffre = re.compile(r"([^;]+)\=([^;]+)")
        self._versionRe = re.compile(r"\.\d+$")
        self._features = set(features)
        self.genes = {}

    def _parser(self) -> Callable:
        """
        Decide which parser to return
        """
        with open(self.gff3, "rb") as _rb:
            if _rb.read(2) == b"\x1f\x8b":
                return partial(gzip.open, mode="rt")
            return partial(open, mode="r")

    def _size(self) -> int:
        """
        Helper function
        Return file size
        """
        with self._reader(self.gff3, mode="rb") as _sh:
            _sh.seek(0, 2)
            size: int = _sh.tell()
        return size

    def parse_gene_info(
        self,
        gene_id: str = "gene_id",
        gene_name: str = "gene_name",
        gene_type="gene_type",
    ) -> Dict[str, List[Gene]]:
        """
        Return gene data
        """
        with self._reader(self.gff3) as fh:
            for l in fh:
                if l[0] == "#":
                    continue
                ldat: List[str] = l.strip().split("\t")
                if (len(ldat) < 9) or (ldat[2] not in self._features):
                    continue
                attribs: Dict[str, str] = dict(re.findall(self._gffre, ldat[8]))
                if gene_id not in attribs:
                    continue
                attribs[gene_id] = attribs[gene_id].split(".")[0]
                try:
                    begin = int(ldat[3]) - 1
                except ValueError:
                    continue
                end = int(ldat[4])
                try:
                    name = attribs[gene_name]
                except KeyError:
                    name = attribs[gene_id]
                try:
                    g_type = attribs[gene_type]
                except KeyError:
                    g_type = "Unknown"
                gene_dat = Gene(
                    attribs[gene_id], name, g_type, ldat[0], begin, end, ldat[6]
                )
                try:
                    self.genes[ldat[0]].append(gene_dat)
                except KeyError:
                    self.genes[ldat[0]] = [gene_dat]
        if len(self.genes) == 0:
            raise RuntimeError(
                f"Cannot parse features: {', '.join(self._features)} from {self.gff3}! Check your inputs!"
            )
        return self.genes
