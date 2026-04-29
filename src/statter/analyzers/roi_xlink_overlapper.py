from collections import namedtuple
from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp
from typing import Callable, Optional

import pandas as pd

from statter.parsers.bed_ROI_parser import RegionReader
from statter.parsers.crosslink_parser import SitesCounter
from statter.parsers.csv_meta_parser import MetaReader


class RegionXlinkOverlapFinder:

    def __init__(
        self,
        metadata: str | Path,
        region: str | Path,
        l: int = 100,
        r: int = 100,
        unstranded: bool = False,
        most_5prime: bool = False,
        norm_method: Optional[str] = "cpm",
        tmpdir: Optional[str | Path] = None,
    ) -> None:
        self.metadata = metadata
        self.region = region
        self.l = l
        self.r = r
        self.unstranded = unstranded
        self.most_5prime = most_5prime
        self.norm_method = norm_method
        self.tmpdir = tmpdir

    def __enter__(self) -> "RegionXlinkOverlapFinder":
        # create temp dir
        if self.tmpdir is None:
            self._tmp = mkdtemp()
        else:
            self._tmp = Path(self.tmpdir).resolve()
            self._tmp.mkdir(exist_ok=True, parents=True)
            self._tmp = mkdtemp(dir=self._tmp)
        # read metadata
        meta_reader = MetaReader(self.metadata)
        self._meta_df = meta_reader.read_meta()
        # all sample names
        self._sample_names: set[str] = set(self._meta_df["sample"].unique())
        # read and index regions
        region_indexer = RegionReader(
            self.region,
            l=self.l,
            r=self.r,
            unstranded=self.unstranded,
            most_5prime=self.most_5prime,
        )
        self._region_df = region_indexer.regions
        self._region_max_len = region_indexer.max_length
        # read and aggregate sites
        sites_to_bed = SitesCounter(
            self._meta_df, norm_method=self.norm_method, tempfolder=self._tmp
        )
        self._sites = sites_to_bed.count_crosslinks()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        # remove temp dir
        rmtree(self._tmp)

    def find_overlaps(self):
        dfs: list[pd.DataFrame] = []
        for region in self._region_df.itertuples(index=False):
            filters = [
                ("chrom", "==", region.chrom),
                ("strand", "==", region.strand),
                ("end", ">", region.start_extend),
                ("end", "<=", region.end_extend),
            ]
            counts: pd.DataFrame = pd.read_parquet(
                self._sites, engine="pyarrow", filters=filters
            )
            if counts.shape[0] == 0:
                continue
            counts = counts.pivot_table(
                index="sample", columns="end", values="count", fill_value=0.0
            )
            all_cols: set[int] = set(
                range(region.start_extend + 1, region.end_extend + 1, 1)  # type: ignore
            )
            # fill missing positions with 0
            missing_cols = sorted(all_cols - set(counts.columns))
            if len(missing_cols) > 0:
                counts = pd.concat(
                    [
                        counts,
                        pd.DataFrame(0.0, index=counts.index, columns=missing_cols),
                    ],
                    axis=1,
                )[sorted(all_cols)]
            # fill missing samples with 0
            missing_samples = self._sample_names - set(counts.index)
            if len(missing_samples) > 0:
                counts = pd.concat(
                    [
                        counts,
                        pd.DataFrame(
                            0.0, index=sorted(missing_samples), columns=counts.columns
                        ),
                    ],
                    axis=0,
                ).sort_index()
                dfs.append(self._format_count_df(counts, region))
        if len(dfs) == 0:
            raise ValueError("No overlaps found between regions and crosslink sites!")
        return pd.concat(dfs, axis=0).reset_index(drop=True)

    def _format_count_df(self, counts: pd.DataFrame, region: namedtuple):
        strand_order_fn = self._order_df_fn(region.strand)
        five_prime: pd.DataFrame = strand_order_fn(
            counts.loc[:, region.end + 1 :]
            if region.strand == "-"
            else counts.loc[:, : region.start]
        )
        five_prime.columns = list(
            range(-(five_prime.shape[1] - 1), 1)
        )  # 5' postitions relative to ROI
        roi: pd.DataFrame = strand_order_fn(
            counts.loc[:, region.start + 1 : region.end]
        )  # bed index
        # fill missing positions with 0s
        if roi.shape[1] < self._region_max_len:
            # pad with 0 to make all regions the same length
            diff = self._region_max_len - roi.shape[1]
            cols = [roi.columns[-1]] * diff
            roi = pd.concat(
                [roi, pd.DataFrame(0.0, index=roi.index, columns=cols)], axis=1
            )
        roi.columns = list(range(1, roi.shape[1] + 1))
        three_prime: pd.DataFrame = strand_order_fn(
            counts.loc[:, : region.start]
            if region.strand == "-"
            else counts.loc[:, region.end + 1 :]
        )
        three_prime.columns = list(
            range(roi.shape[1] + 1, roi.shape[1] + three_prime.shape[1] + 1)
        )  # 3' positions relative to ROI
        return pd.concat([five_prime, roi, three_prime], axis=1)

    def _order_df_fn(self, strand: str) -> Callable[[pd.DataFrame], pd.DataFrame]:

        def _plus_strand_order(df: pd.DataFrame) -> pd.DataFrame:
            return df

        def _minus_strand_order(df: pd.DataFrame) -> pd.DataFrame:
            return df.iloc[:, ::-1]

        if strand == "-":
            return _minus_strand_order
        else:
            return _plus_strand_order
