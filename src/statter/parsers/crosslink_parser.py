import logging
import tempfile
from pathlib import Path
from typing import Callable, Optional

import pandas as pd

logger = logging.getLogger(__name__)


class SitesCounter:
    """
    convert Shoji/htseq-clip sites data to BED where the score column has the cpm normalized count value
    """

    def __init__(
        self,
        metadf: pd.DataFrame,
        norm_method: Optional[str] = None,
        tempfolder: Optional[str | Path] = None,
    ) -> None:

        self.metadata = metadf
        self.norm = norm_method
        if tempfolder is None:
            self._tmpdir = Path(tempfile.gettempdir())
        else:
            self._tmpdir = Path(tempfolder)
            self._tmpdir.mkdir(exist_ok=True, parents=True)

        if self.norm is not None:
            self.norm = self.norm.lower()
        # for bed file
        self._cols: list[str] = ["chrom", "start", "end", "name", "score", "strand"]
        self._count_cols: list[str] = ["chrom", "end", "strand"]
        self._partition_cols: list[str] = ["chrom", "strand"]
        self._normalizer = self._get_normalizer()

    def count_crosslinks(self) -> Path:
        pp_df: Path = self._tmpdir / next(tempfile._get_candidate_names())  # type: ignore
        for dat in self.metadata.itertuples():
            logging.debug(
                "Processing file {} for sample {} in group {}".format(
                    dat.file, dat.sample, dat.group
                )
            )
            xbed = self._xlink_reader(dat.file, dat.sample, dat.group)
            xbed.to_parquet(
                pp_df,
                engine="pyarrow",
                partition_cols=self._partition_cols,
                index=False,
            )
        return pp_df

    def _get_normalizer(self) -> Callable[[pd.DataFrame], pd.DataFrame]:
        """_get_normalizer Helper function
        Return appropriate function to normalize data

        Returns:
            Callable[[pd.DataFrame], pd.DataFrame]: function to normalize data
        """

        def donothing(xbed: pd.DataFrame) -> pd.DataFrame:
            return xbed

        def cpm(xbed: pd.DataFrame) -> pd.DataFrame:
            xbed["count"] = xbed["count"] * 1e6 / xbed["count"].sum()
            return xbed

        if self.norm is None:
            return donothing
        elif self.norm == "cpm":
            return cpm
        else:
            raise NotImplementedError(
                "Given normalization method {} is not implemented".format(self.norm)
            )

    def _xlink_reader(self, fname: str, sample: str, group: str) -> pd.DataFrame:
        """_xlink_reader Read and aggregate crosslink sites from bed file

        Args:
            fname: (str) path to bed file containing crosslink sites
            sample: (str) name of the sample
            group: (str) group to which the sample belongs

        Returns:
            pd.DataFrame: DataFrame containing the aggregated and optionally normalized crosslink sites
        """
        xbed: pd.DataFrame = (
            pd.read_csv(fname, sep="\t", header=None, names=self._cols)
            .groupby(self._count_cols)
            .size()
            .reset_index(name="count")
        )
        xbed["sample"] = sample
        xbed["group"] = group
        return self._normalizer(xbed)

    # def aggregate(self):
    #     bedmap = {}
    #     for fname, dat in self.filedict.items():
    #         freader, mode = self._get_parser(fname)
    #         countdict = {}
    #         with freader(fname, mode=mode) as fh:
    #             logger.info("Aggregating crosslink counts from {}".format(fname))
    #             for f in fh:
    #                 f = f.strip().split("\t")
    #                 try:
    #                     countdict[(f[0], f[1], f[2], f[-1])] += 1
    #                 except KeyError:
    #                     countdict[(f[0], f[1], f[2], f[-1])] = 1
    #         out_name = re.sub(r"\s{1,}", "_", dat["name"])  # avoid spaces
    #         sout = str(self.tempdir / "{}.bed".format(out_name))
    #         count_sum = np.sum(list(countdict.values()))
    #         bedmap[sout] = dat
    #         with open(sout, "w") as oh:
    #             logger.debug("Writing {} to {}".format(fname, sout))
    #             for k, v in countdict.items():
    #                 oh.write(
    #                     "{}\t{}\t{}\t{}\t{}\t{}\n".format(
    #                         k[0],
    #                         k[1],
    #                         k[2],
    #                         dat["name"],
    #                         self._normalizer(v, count_sum),  # type: ignore
    #                         k[3],
    #                     )
    #                 )
    #     return bedmap
