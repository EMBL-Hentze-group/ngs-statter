import gzip
import logging
import re
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


class SitesToBed:
    """
    convert htseq-clip sites data to BED where the score column has the cpm normalized count value
    """

    def __init__(self, tempfolder, filedict, norm_method=None) -> None:
        self.tempdir = Path(tempfolder)
        self.filedict = filedict
        if not self.tempdir.exists():
            logger.debug("Creating temp dir {}".format(self.tempdir))
            self.tempdir.mkdir()
        self.norm = norm_method
        if self.norm is not None:
            self.norm = self.norm.lower()
        self._norm_methods = set(["cpm"])
        self._normalizer = self._get_normalizer()

    @staticmethod
    def _get_parser(fname):
        with open(fname, "rb") as _rb:
            if _rb.read(2) == b"\x1f\x8b":
                return (gzip.open, "rt")
            else:
                return (open, "r")

    def _get_normalizer(self):

        def donothing(c, csum):
            return c

        def cpm(c, csum):
            return float(c / csum) * 1e6

        if self.norm is None:
            return donothing
        else:
            if self.norm not in self._norm_methods:
                raise NotImplementedError(
                    "Given normalization method {} is not implemented".format(self.norm)
                )
            elif self.norm == "cpm":
                return cpm

    def aggregate(self):
        bedmap = {}
        for fname, dat in self.filedict.items():
            freader, mode = self._get_parser(fname)
            countdict = {}
            with freader(fname, mode=mode) as fh:
                logger.info("Aggregating crosslink counts from {}".format(fname))
                for f in fh:
                    f = f.strip().split("\t")
                    try:
                        countdict[(f[0], f[1], f[2], f[-1])] += 1
                    except KeyError:
                        countdict[(f[0], f[1], f[2], f[-1])] = 1
            out_name = re.sub(r"\s{1,}", "_", dat["name"])  # avoid spaces
            sout = str(self.tempdir / "{}.bed".format(out_name))
            count_sum = np.sum(list(countdict.values()))
            bedmap[sout] = dat
            with open(sout, "w") as oh:
                logger.debug("Writing {} to {}".format(fname, sout))
                for k, v in countdict.items():
                    oh.write(
                        "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            k[0],
                            k[1],
                            k[2],
                            dat["name"],
                            self._normalizer(v, count_sum),  # type: ignore
                            k[3],
                        )
                    )
        return bedmap
