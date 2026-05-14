import logging
from pathlib import Path
from typing import Generator

import polars as pl

logger = logging.getLogger(__name__)


class RegionReader:

    def __init__(
        self,
        bed: str | Path,
        l: int = 100,
        r: int = 100,
        unstranded: bool = False,
        most_5prime: bool = False,
    ) -> None:
        self._schema = {
            "chrom": pl.String,
            "start": pl.UInt32,
            "end": pl.UInt32,
            "name": pl.String,
            "score": pl.Float32,
            "strand": pl.String,
        }
        self._group_cols = ["chrom", "strand"]
        self._lf = pl.scan_csv(
            bed,
            separator="\t",
            has_header=False,
            schema=self._schema,
            comment_prefix="#",
        )
        self._check_nulls()
        self._slop(l, r, unstranded)
        if most_5prime:
            self._pick_5prime_most()
        # renaming to avoid confusion with the "end" column in crosslink
        # casting to int64 to avoid overflow issues in downstream calculations
        self._lf = self._lf.rename({"end": "stop"}).cast(
            {"start": pl.Int64, "stop": pl.Int64}
        )
        self._max_length: int = (
            self._lf.select((pl.col("stop") - pl.col("start")).max()).collect().item()
        )

    def _check_nulls(self) -> None:
        """
        check if any of the required columns contain null values.
        Raises:
            ValueError: If any required column contains null values.
        """
        non_nulls = (
            self._lf.select(pl.all().null_count() == 0).collect().to_dicts()[0].values()
        )
        if not all(non_nulls):
            raise ValueError("BED file contains null values in required columns.")

    def _slop(self, l: int, r: int, unstranded: bool) -> None:
        """_slop Helper function
        Extend regions in 5' and 3' directions by given lengths l and r.
        Args:
            l: (int) number of bases to extend at the 5' end
            r: (int) number of bases to extend at the 3' end
            unstranded: (bool) whether to ignore strand information
        """
        if unstranded:
            self._lf = self._lf.with_columns(
                [
                    (pl.col("start") - l).clip(0).alias("start_extend"),
                    (pl.col("end") + r).alias("end_extend"),
                ]
            )
        else:
            self._lf = self._lf.with_columns(
                [
                    pl.when(pl.col("strand") == "-")
                    .then((pl.col("start") - r).clip(0))
                    .otherwise((pl.col("start") - l).clip(0))
                    .alias("start_extend"),
                    pl.when(pl.col("strand") == "-")
                    .then(pl.col("end") + l)
                    .otherwise(pl.col("end") + r)
                    .alias("end_extend"),
                ]
            )

    def _pick_5prime_most(self):
        """_pick_5prime_most Helper function

        Pick the most 5' regions out of the set of overlapping regions after extension
        """
        sort_cols: list[str] = self._group_cols + ["start_extend", "end_extend"]
        self._lf = (
            self._lf.sort(sort_cols)
            .with_columns(
                [
                    pl.when(pl.col("strand") == "-")
                    .then(
                        pl.col("start_extend")
                        .cum_max()
                        .shift(-1)
                        .over(self._group_cols)
                    )
                    .otherwise(
                        pl.col("end_extend").cum_max().shift(1).over(self._group_cols)
                    )
                    .alias("prev_boundary")
                ]
            )
            .filter(
                (pl.col("prev_boundary").is_null())
                | (
                    (pl.col("strand") == "+")
                    & (pl.col("start_extend") >= pl.col("prev_boundary"))
                )
                | (
                    (pl.col("strand") == "-")
                    & (pl.col("end_extend") <= pl.col("prev_boundary"))
                )
            )
            .drop("prev_boundary")
        )

    @property
    def regions(self) -> pl.LazyFrame:
        return self._lf

    @property
    def max_length(self) -> int:
        """
        return the max. length of the ROIs before extension
        """
        return self._max_length

    def filter_queries(
        self, min_gap: int = 50000, batch_size: int = 1000
    ) -> Generator[pl.DataFrame, None, None]:
        """filter_queries
        This function yields batches of regions of interest (ROIs) after extension, where the gap between consecutive ROIs on the same chromosome and strand is at least `min_gap`. The dataframe is yielded in batches of size `batch_size` to manage memory for large datasets.
        This is used in downstream steps to filter crosslink sites that fall within these ROIs.
        Args:
            min_gap: (int) minimum gap between consecutive ROIs. Defaults to 50000.
            batch_size: (int) number of ROIs to yield per batch. Defaults to 1000.

        Yields:
            pl.DataFrame: A batch of consecutive ROIs.
        """
        query: pl.LazyFrame = (
            self._lf.sort(["chrom", "strand", "end_extend"])
            .with_columns(gap=pl.col("end_extend").diff().over(["chrom", "strand"]))
            .with_columns(
                end_extend=pl.when(pl.col("gap").is_null())
                .then(pl.col("start_extend"))
                .otherwise(pl.col("end_extend"))
            )
            .filter(pl.col("gap").is_null() | (pl.col("gap") >= min_gap))
            .with_columns(
                prev_extend=pl.col("end_extend").shift(1).over(["chrom", "strand"])
            )
            .filter(~pl.col("gap").is_null())
            .select(["chrom", "strand", "prev_extend", "end_extend"])
        )
        for batch in query.collect_batches(chunk_size=batch_size):
            yield batch
