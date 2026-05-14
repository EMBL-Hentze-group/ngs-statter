import logging
import tempfile
from pathlib import Path
from typing import Callable, Optional

import polars as pl

logger = logging.getLogger(__name__)


class SitesCounter:

    def __init__(
        self,
        metadf: pl.DataFrame,
        norm_method: Optional[str] = None,
        tempfolder: Optional[str | Path] = None,
    ) -> None:
        """__init__

        Args:
            metadf: (pl.DataFrame) metadata dataframe containing columns 'file', 'sample', and 'group'
            norm_method: (Optional[str]) normalization method to apply. Defaults to None.
            tempfolder: (Optional[str | Path]) path to temporary folder. Defaults to None.
        """
        self.metadata = metadf
        self._norm_fn = self._normalizer(norm_method)
        if tempfolder is None:
            self._tmpdir = Path(tempfile.gettempdir())
        else:
            self._tmpdir = Path(tempfolder)
            self._tmpdir.mkdir(exist_ok=True, parents=True)
        # for bed file
        self._bed_schema: dict[str, pl.DataType] = {  # type: ignore
            "chrom": pl.String,
            "start": pl.UInt32,
            "end": pl.Int64,  # to avoid integer overflow during join in the next steps
            "name": pl.String,
            "score": pl.Float32,
            "strand": pl.String,
        }
        self._group_by_cols: list[str] = ["chrom", "end", "strand"]
        self._partition_cols: list[str] = ["chrom", "strand"]
        # for filtering
        self._groups: list[str] = (
            self.metadata.get_column("group").unique().sort().to_list()
        )
        self._samples: list[str] = (
            self.metadata.get_column("sample").unique().sort().to_list()
        )

    def _normalizer(
        self, norm: Optional[str]
    ) -> Callable[[pl.LazyFrame, int], pl.LazyFrame]:
        """_normalizer Helper function
        Return a normalization function based on the specified normalization method.
        Args:
            norm: (Optional[str]) normalization method to apply.

        Returns:
            Callable[[pl.LazyFrame, int], pl.LazyFrame]: A function that takes a LazyFrame and library size, and returns a normalized LazyFrame.
        """

        def nothing(df: pl.LazyFrame, lib_size: int) -> pl.LazyFrame:
            return df

        def cpm(df: pl.LazyFrame, lib_size: int) -> pl.LazyFrame:
            return df.with_columns((pl.col("count") / lib_size * 1e6).alias("count"))

        if norm is None:
            return nothing
        norm = norm.lower()
        if norm == "cpm":
            return cpm
        else:
            raise NotImplementedError(
                f"Normalization method {norm} is not implemented."
            )

    def count_crosslinks(self) -> Path:
        """count_crosslinks
        Aggregate crosslink sites per sample, write to a paritioned parquet file and return file path

        Returns:
            Path to the partitioned parquet file
        """
        pp_df: Path = self._tmpdir / next(tempfile._get_candidate_names())  # type: ignore
        pp_df.mkdir(exist_ok=True, parents=True)
        for dat in self.metadata.iter_rows(named=True):
            logging.info(
                f"Processing file {dat['file']} for sample {dat['sample']} in group {dat['group']}"
            )
            self._xlink_parser(dat["file"], dat["sample"], dat["group"], pp_df)
        return pp_df

    def _xlink_parser(
        self, filepath: str, sample: str, group: str, out_folder: Path
    ) -> None:
        """_xlink_parser Helper function
        Parse crosslink sites per sample, aggregate counts, normalize, and write to partitioned parquet directory
        Args:
            filepath: (str) path to the input file
            sample: (str) sample name
            group: (str) group name
            out_folder: (Path) path to the output folder
        """
        xlinks: pl.LazyFrame = (
            pl.scan_csv(
                filepath,
                has_header=False,
                separator="\t",
                comment_prefix="#",
                schema=self._bed_schema,
            )
            .group_by(self._group_by_cols)
            .agg(pl.len().alias("count").cast(pl.Float32))
            .with_columns(
                [pl.lit(sample).alias("sample"), pl.lit(group).alias("group")]
            )
        )

        lib_size: int = xlinks.select(pl.sum("count")).collect().item()
        xlinks = self._norm_fn(xlinks, lib_size)
        xlinks.sink_parquet(
            pl.PartitionBy(
                base_path=out_folder,
                file_path_provider=_temp_name_provider,
                key=self._partition_cols,
                include_key=False,
            ),
            mkdir=True,
        )


def _temp_name_provider(args) -> str:
    """_temp_name_provider Helper function
    Generate a temporary file name for partitioned parquet files based on the partition keys.
    see
    https://docs.pola.rs/api/python/stable/reference/api/polars.LazyFrame.sink_parquet.html and
    https://docs.pola.rs/api/python/stable/reference/api/polars.PartitionBy.html#polars.PartitionBy
    Args:
        args: Arguments containing partition keys.

    Returns:
        str: Path to the temporary parquet file.
    """
    fps: list[str] = []
    df = args.partition_keys
    for col in df.columns:
        fps.append(f"{col}={df[col][0]}")
    path = Path().joinpath(*fps) / f"{next(tempfile._get_candidate_names())}.parquet"  # type: ignore
    return str(path)
