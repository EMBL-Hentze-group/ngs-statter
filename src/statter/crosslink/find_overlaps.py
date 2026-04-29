from pathlib import Path
from shutil import rmtree
from tempfile import mkdtemp
from typing import Optional

import pandas as pd

from statter.parsers.crosslink.meta_parser import MetaReader
from statter.parsers.crosslink.sites_to_bed import SitesCounter
from statter.parsers.crosslink.structure_indexer import RegionReader


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
        self._sample_names: set[str] = set(self._meta_df["sample_name"].unique())
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
            missing_cols = all_cols - set(counts.columns)
            if len(missing_cols) > 0:
                counts = pd.concat(
                    [
                        counts,
                        pd.DataFrame(0.0, index=counts.index, columns=missing_cols),
                    ],
                    axis=1,
                )[sorted(all_cols)]
            missing_samples = self._sample_names - set(counts.index)
            if len(missing_samples) > 0:
                counts = pd.concat(
                    [
                        counts,
                        pd.DataFrame(
                            0.0, index=missing_samples, columns=counts.columns
                        ),
                    ],
                    axis=0,
                ).sort_index()
