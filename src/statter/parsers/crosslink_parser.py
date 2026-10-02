import logging
from collections.abc import Callable
from pathlib import Path

import polars as pl

from statter.parsers.base_data_parser import BaseDataFrameProcessor

logger = logging.getLogger(__name__)


class SitesCounter(BaseDataFrameProcessor):
    def __init__(
        self,
        metadf: pl.DataFrame,
        norm_method: str | None,
        tempfolder: str | Path | None,
    ) -> None:
        super().__init__(metadf, tempfolder=tempfolder)
        self._norm_fn = self._normalizer(norm_method)

    @property
    def _schema(self) -> dict[str, pl.DataType]:
        return {
            "chrom": pl.String,
            "start": pl.UInt32,
            "end": pl.Int64,
            "name": pl.String,
            "score": pl.Float32,
            "strand": pl.String,
        }  # type: ignore

    @property
    def _group_by_cols(self) -> list[str]:
        return ["chrom", "end", "strand"]

    @property
    def _partition_cols(self) -> list[str]:
        return ["chrom", "strand"]

    @property
    def _rename_cols(self) -> dict[str, str] | None:
        return None

    def _normalizer(
        self, norm: str | None
    ) -> Callable[[pl.LazyFrame, int], pl.LazyFrame]:
        def nothing(df: pl.LazyFrame, lib_size: int) -> pl.LazyFrame:
            return df

        def cpm(df: pl.LazyFrame, lib_size: int) -> pl.LazyFrame:
            return df.with_columns((pl.col("count") / lib_size * 1e6).alias("count"))

        if norm is None:
            return nothing
        norm = norm.lower()
        if norm == "cpm":
            return cpm
        raise NotImplementedError(f"Normalization method {norm} is not implemented.")

    def _parse_and_transform(
        self, filepath: str, sample: str, group: str
    ) -> pl.LazyFrame:
        xlinks = (
            pl.scan_csv(
                filepath,
                has_header=False,
                separator="\t",
                comment_prefix="#",
                schema=self._schema,
            )
            .group_by(self._group_by_cols)
            .agg(pl.len().alias("count").cast(pl.Float32))
            .with_columns(
                [
                    pl.lit(sample).alias("sample"),
                    pl.lit(group).alias("group"),
                ]
            )
        )
        if self._rename_cols is not None:
            xlinks = xlinks.rename(self._rename_cols)
        lib_size = xlinks.select(pl.sum("count")).collect().item()
        return self._norm_fn(xlinks, lib_size)
