import logging
from pathlib import Path

import polars as pl

from statter.analyzers.roi_xlink_overlapper import RegionSignalOverlapFinder
from statter.parsers.allele_depth_parser import AlleleDepthParser
from statter.parsers.crosslink_parser import SitesCounter
from statter.parsers.csv_meta_parser import MetaReader

logger = logging.getLogger(__name__)


class RoiOverlapWrapper:
    def __init__(
        self,
        meta_data: str | Path,
        region: str | Path,
        l: int,
        r: int,
        unstranded: bool,
        most_5prime: bool,
        smoothing_window: int,
        out_file: str,
        tmp_dir: str | Path | None,
    ):
        self.meta_data = meta_data
        self.region = region
        self.l = l
        self.r = r
        self.unstranded = unstranded
        self.most_5prime = most_5prime
        self.smoothing_window = smoothing_window
        self.out_file = out_file
        self.tmp_dir = tmp_dir
        # params to be set by metadata parser
        self._meta_df: pl.DataFrame  # lazy frame containing metadata information
        self._group_col_df: (
            pl.LazyFrame
        )  # lazy frame containing group color information
        # params to be set by signal processing class
        self._signal_dir: Path  # directory containing signal files
        self._partition_cols: list[
            str
        ]  # list of columns used for partitioning the signal data
        #  parse metadata
        self._parse_meta()
        # max length for the region
        self._region_max_len: int  # maximum length of the region

    @property
    def region_max_len(self) -> int:
        return self._region_max_len

    def _parse_meta(self):
        """_parse_meta Helper function
        Parses the metadata file and initializes the metadata DataFrame and group color DataFrame.
        """
        meta_reader = MetaReader(self.meta_data)
        self._meta_df = meta_reader.read_meta()
        self._group_col_df = meta_reader.per_group_color_lf().lazy()

    def _overlapper(self):
        with RegionSignalOverlapFinder(
            signal_dir=self._signal_dir,
            partition_cols=self._partition_cols,
            group_color=self._group_col_df,
            region=self.region,
            l=self.l,
            r=self.r,
            unstranded=self.unstranded,
            most_5prime=self.most_5prime,
            smoothing_window=self.smoothing_window,
            tmpdir=self._signal_dir.parent,
        ) as ros:
            ros.find_overlaps(self.out_file)
            self._region_max_len = ros.region_max_len

    def crosslink_overlap(self, norm_method: str | None):
        xov = SitesCounter(
            metadf=self._meta_df, norm_method=norm_method, tempfolder=self.tmp_dir
        )
        self._signal_dir = xov.to_signal()
        print(f"Signal directory: {self._signal_dir}")

        self._partition_cols = xov.partition_cols
        self._overlapper()

    def allele_depth_overlap(
        self,
        genome_fa: str | Path | None,
        variant_type: str | None,
        dp: int = 0,
        ncores: int = 1,
    ):
        adov = AlleleDepthParser(
            metadf=self._meta_df,
            tempfolder=self.tmp_dir,
            genome_fa=genome_fa,
            variant_type=variant_type,
            dp=dp,
            ncores=ncores,
        )
        self._signal_dir = adov.to_signal()
        self._partition_cols = adov.partition_cols
        self._overlapper()
