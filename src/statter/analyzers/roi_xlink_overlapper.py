import logging
import tempfile
from pathlib import Path
from shutil import rmtree
from typing import Optional, Self

import numpy as np
import pandas as pd
import polars as pl
from scipy.ndimage import convolve1d
from scipy.sparse import coo_matrix

from statter.parsers.bed_ROI_parser import RegionReader

logger = logging.getLogger(__name__)


class RegionSignalOverlapFinder:

    def __init__(
        self,
        signal_dir: str | Path,
        partition_cols: list[str],
        group_color: pl.LazyFrame,
        region: str | Path,
        # mode: str,
        l: int = 100,
        r: int = 100,
        unstranded: bool = False,
        most_5prime: bool = False,
        smoothing_window: int = 0,
        # norm_method: str | None = "cpm",
        tmpdir: str | Path | None = None,
    ) -> None:
        """__init__ _summary_

        Args:
            signal_dir: (str or Path) path containing partitioned signal files (.parquet format)
            partition_cols: (list[str]) columns used to partition the signal data for efficient processing.
            group_color: (pl.LazyFrame) lazy frame containing group-color mapping. Has two columns: 'group' and 'color'.
            region: (str or Path) path to region file (bed)
            l: (int) left extension of the region. Defaults to 100.
            r: (int) right extension of the region. Defaults to 100.
            unstranded: (bool) whether to ignore strand information. Defaults to False.
            most_5prime: (bool) whether to consider the most 5' regions in overlaps. Defaults to False.
            smoothing_window: (int) window size for smoothing. Defaults to 0.
            tmpdir: (Optional[str or Path]) path to temporary directory. Defaults to None.
        """
        self.signal_dir = signal_dir
        self.region = region
        self.group_color = group_color
        self.partition_cols = partition_cols
        self.l = l
        self.r = r
        # self.mode = mode
        self.unstranded = unstranded
        self.most_5prime = most_5prime
        self.smoothing_window = smoothing_window
        # currently only "cpm" normalization method is supported, for crosslink plot visualization
        # self.norm_method = norm_method
        self.tmpdir = self._init_tmpdir(tmpdir)
        # read metadata
        # self._meta_reader: MetaReader = MetaReader(self.metadata)
        # self._meta_df: pl.DataFrame = self._meta_reader.read_meta()
        # signals
        # self._signal: AlleleDepthParser | SitesCounter | None = None
        # use these columns to partition the data for efficient processing
        # signal parser and partitioning columns are initialized in _get_signal()
        # self._get_signal()

    def _init_tmpdir(self, tempfolder: str | Path | None) -> Path:
        """Initialize the temporary directory.
        Args:
            tempfolder: (str or Path or None) path to the temporary directory. If None, the system temporary directory is used.

        Returns:
            Path to the initialized temporary directory.
        """
        if tempfolder is None:
            return Path(tempfile.gettempdir())
        tmpdir = Path(tempfolder)
        tmpdir.mkdir(exist_ok=True, parents=True)
        return tmpdir

    # def _get_signal(self) -> None:
    #     if self.mode == "allele_depth":
    #         self._signal = AlleleDepthParser(self._meta_df, self.tmpdir)
    #         self._by_cols = ["chrom"]
    #     elif self.mode == "crosslink":
    #         self._signal = SitesCounter(self._meta_df, self.norm_method, self.tmpdir)
    #         self._by_cols = ["chrom", "strand"]
    #     else:
    #         raise NotImplementedError(f"Unsupported mode: {self.mode}")

    def __enter__(self) -> Self:
        # read and index regions
        region_indexer = RegionReader(
            self.region,
            l=self.l,
            r=self.r,
            unstranded=self.unstranded,
            most_5prime=self.most_5prime,
        )
        self._region_max_len: int = region_indexer.max_length
        regions: pl.LazyFrame = region_indexer.regions
        # read and aggregate crosslink sites
        # if self._signal is None:
        #     raise ValueError("Signal parser is not initialized.")
        sites: pl.LazyFrame = pl.scan_parquet(self.signal_dir, hive_partitioning=True)
        # define enums for chromosome and strand to optimize joins
        chroms = pl.Enum(
            pl.concat(
                [
                    regions.select("chrom").unique(),
                    sites.select("chrom").unique(),
                ]
            )
            .unique()
            .sort("chrom")
            .collect()
            .to_series()
        )
        strands = pl.Enum(
            pl.concat(
                [
                    regions.select("strand").unique(),
                    sites.select("strand").unique(),
                ]
            )
            .unique()
            .sort("strand")
            .collect()
            .to_series()
        )
        # cast chromosome and strand columns to enums
        regions = regions.with_columns(
            pl.col("chrom").cast(chroms),
            pl.col("strand").cast(strands),
        )
        sites = sites.with_columns(
            pl.col("chrom").cast(chroms),
            pl.col("strand").cast(strands),
        )
        # sink regions and sites to parquet files for efficient querying
        regions: pl.LazyFrame = pl.scan_parquet(self._sink_regions(regions))
        sites: pl.LazyFrame = pl.scan_parquet(self._sink_sites(sites))
        # query to find overlaps
        self._query = (
            sites.join_asof(
                regions,
                left_on="end",
                right_on="start_extend",
                by=self.partition_cols,
                strategy="backward",
                tolerance=(self.l + self.r + self._region_max_len),
                allow_exact_matches=False,
            )
            .filter(pl.col("end") <= pl.col("end_extend"))
            .with_columns(
                rlen=pl.col("end") - pl.col("start"),
                len_diff=(
                    self._region_max_len - (pl.col("stop") - pl.col("start"))
                ).clip(0),
                rel_pos=pl.when(pl.col("strand") == "-")
                .then(pl.col("stop") - pl.col("end"))
                .otherwise(pl.col("end") - pl.col("start")),
            )
            .with_columns(
                rel_pos=pl.when(pl.col("rel_pos") > pl.col("rlen"))
                .then(pl.col("rel_pos") + pl.col("len_diff"))
                .otherwise(pl.col("rel_pos"))
            )
            .drop(["rlen", "len_diff"])
        )
        # define required columns for counting overlaps (positions relative to region)
        self._req_cols: list[str] = list(
            map(str, range(-self.l + 1, self.r + self._region_max_len + 1, 1))
        )
        self._col_set = set(self._req_cols)
        _index_cols: list[str] = [
            "group",
            "sample",
        ]
        # columns containing position details
        _pos_cols: list[str] = [
            "chrom",
            "start",
            "stop",
            "strand",
        ]
        self._index_cols: list[str] = _pos_cols + _index_cols
        return self

    def _get_tmp_path(self, suffix: str | None) -> Path:

        if suffix is None:
            return self.tmpdir / f"{next(tempfile._get_candidate_names())}.parquet"  # type: ignore

        return self.tmpdir / f"{next(tempfile._get_candidate_names())}_{suffix}.parquet"  # type: ignore

    def _sink_regions(self, regions: pl.LazyFrame) -> Path:
        """_sink_regions Helper function
        write regions to a parquet file and return file path

        Args:
            regions: (pl.LazyFrame) regions to be written to parquet file

        Returns:
            Path: file path of the written parquet file
        """
        file_path: Path = self._get_tmp_path("regions")
        regions.drop(["name", "score"]).sort(
            ["chrom", "strand", "start_extend", "end_extend"]
        ).sink_parquet(file_path)
        return file_path

    def _sink_sites(self, sites: pl.LazyFrame) -> Path:
        """_sink_sites Helper function
        write crosslink sites to a parquet file and return file path

        Args:
            sites: (pl.LazyFrame) crosslink sites to be written to parquet file

        Returns:
            Path: file path of the written parquet file
        """
        file_path: Path = self._get_tmp_path("sites")
        sites.sort(["chrom", "strand", "end"]).sink_parquet(file_path)
        return file_path

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        rmtree(self.tmpdir)

    @property
    def region_max_len(self) -> int:
        return self._region_max_len

    def find_overlaps(self, out_path: str) -> None:
        """find_overlaps
        Find overlaps between the given region(s) of interest and crosslink sites
        write the data to specified output path in parquet format with columns for group, sample and relative position to the region.
        Args:
            out_path: (str) Path to the output parquet file.
        """

        hits: pl.LazyFrame = self._query.with_columns(
            pl.struct(self._index_cols).rank("dense").alias("row_id") - 1,
            pl.col("rel_pos").rank("dense").alias("col_id") - 1,
        )
        index_path, count_path = self._process_hits(hits)
        # merge index and count files
        if Path(out_path).exists():
            logger.warning(
                f"Output file {out_path} already exists and will be overwritten."
            )

        pl.concat(
            [pl.scan_parquet(index_path), pl.scan_parquet(count_path)], how="horizontal"
        ).sink_parquet(out_path)

    def _process_hits(self, hits: pl.LazyFrame) -> tuple[Path, Path]:
        """_process_hits Helper function
        Process hits to counts and index files
        Args:
            hits: (pl.LazyFrame) hits to be processed

        Returns:
            _tuple[Path, Path]: file paths
                            : first file contains index columns (group, sample, chrom, start, stop, strand)
                            : second file contains counts
        """
        # convert hits to sparse matrix format and apply smoothing
        counts = coo_matrix(
            (
                hits.select("count").collect().to_numpy().flatten(),
                (
                    hits.select("row_id").collect().to_numpy().flatten(),
                    hits.select("col_id").collect().to_numpy().flatten(),
                ),
            )
        )
        sparse_mat: pd.DataFrame = pd.DataFrame(self._smoother(counts.todense(), self.smoothing_window))  # type: ignore
        sparse_mat.columns = (
            hits.select(["col_id", "rel_pos"])
            .unique()
            .sort("col_id")
            .collect()
            .get_column("rel_pos")
            .to_list()
        )
        # save sparse matrix
        sparse_path: Path = self._get_tmp_path("counts")  # type: ignore
        sparse_mat.to_parquet(sparse_path)
        # save index columns (group, sample, chrom, start, stop, strand) and add group colors
        index_path = self._get_tmp_path("index")  # type: ignore
        hits.select(self._index_cols + ["row_id"]).unique().join(
            self.group_color, on="group", how="left"
        ).with_columns(
            pl.col("chrom").cast(pl.String),
            pl.col("strand").cast(pl.String),
        ).sort(
            "row_id"
        ).drop(
            "row_id"
        ).sink_parquet(
            index_path
        )
        return index_path, sparse_path

    # def _get_group_colors(self) -> pl.LazyFrame:
    #     """_get_group_colors Helper function
    #     Get the colors for each group.

    #     Returns:
    #         pl.LazyFrame: DataFrame with columns 'group' and 'color'.
    #     """
    #     groups: list[str] = sorted(self._meta_reader.per_group_colors.keys())
    #     colors: list[str] = [self._meta_reader.per_group_colors[g] for g in groups]
    #     return pl.DataFrame({"group": groups, "color": colors}).lazy()

    def _smoother(self, counts, smoothing_window: int) -> np.ndarray:
        """_smoother Helper function
        Group dataframe by group and sample, then apply smoothing to the values. The smoothing is done by convolving the values with a kernel.
        Args:
            df: (pd.DataFrame) DataFrame to be grouped and smoothed. Must contain 'group' and 'sample' columns.
            smoothing_window: (int) Size of the smoothing window. If 2, a simple average of the current and next value is taken. If greater than 2, a Gaussian kernel is used for smoothing.

        Returns:
            pd.DataFrame: The grouped and smoothed DataFrame.
        """

        if smoothing_window < 2:
            return counts
        kernel = self._kernel(smoothing_window)
        return convolve1d(counts, kernel, axis=1, mode="constant", cval=0.0)

    def _kernel(self, smoothing_window: int) -> np.ndarray:
        """_kernel Helper function
        Generate a smoothing kernel
        Args:
            smoothing_window: (int) Size of the smoothing window.

        Returns:
            np.ndarray: The smoothing kernel.
        """
        if smoothing_window == 2:
            return np.array([0.5, 0.5])
        else:
            k = np.arange(-smoothing_window // 2 + 1, smoothing_window // 2 + 1, 1)
            kernel = np.exp(-(k**2) / (2.0**2))
            return kernel / kernel.sum()
