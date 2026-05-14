import logging
import re
from pathlib import Path
from typing import Optional

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns

logger = logging.getLogger(__name__)


class HeatMapper:

    def __init__(self, df_path: Path | str) -> None:
        self._req_cols: list[str] = [
            "chrom",
            "start",
            "strand",
            "group",
            "sample",
            "color",
        ]
        self._pos_col: str = "position"
        self._count_col: str = "count"
        self._group_cols: list[str] = [
            "chrom",
            "strand",
            "start",
            "group",
            "color",
            "position",
        ]
        self.df_path = df_path
        self._df: pl.LazyFrame = pl.scan_parquet(df_path)
        self._col_sanity_check()
        self._rel_pos: list[str] = self._get_rel_pos_cols()

    def _col_sanity_check(self) -> None:
        missing_cols: set[str] = set(self._req_cols) - set(
            self._df.collect_schema().names()
        )
        if missing_cols:
            raise RuntimeError(
                f"The following columns are missing from {self.df_path}: {' ,'.join(missing_cols)}!"
            )

    def _get_rel_pos_cols(self) -> list[str]:
        rel_pos: list[str] = list(
            filter(
                lambda c: re.match(r"^\-*\d+$", c), self._df.collect_schema().names()
            )
        )
        if len(rel_pos) == 0:
            raise RuntimeError(
                f"{self.df_path} does not contain any relative position columns!"
            )
        return rel_pos

    def _get_cmap(
        self, end_col: str, base_col: str = "#ffffff", name: str = "cmap"
    ) -> mcolors.LinearSegmentedColormap:
        return mcolors.LinearSegmentedColormap.from_list(name, [base_col, end_col])

    def plot(
        self,
        outdir: str,
        width: float = 10,
        height: float = 6,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        roi_length: Optional[int] = None,
        xlabel: str = "Relative position",
        ylabel: str = "Crosslink counts",
    ) -> None:
        out_dir = Path(outdir)
        out_dir.mkdir(parents=True, exist_ok=True)
        df: pl.LazyFrame = (
            self._df.unpivot(
                self._rel_pos,
                index=self._req_cols,
                variable_name=self._pos_col,
                value_name=self._count_col,
            )
            .with_columns(pl.col(self._pos_col).cast(pl.Int32))
            .group_by(self._group_cols)
            .agg(pl.col(self._count_col).mean())
            .sort(self._group_cols)
        )
        if vmin is None:
            vmin = df.select(self._count_col).min().collect().item()
        if vmax is None:
            vmax = df.select(self._count_col).max().collect().item()
        groups: list[str] = (
            df.select("group").unique().collect().get_column("group").to_list()
        )
        for g in groups:
            init: pl.LazyFrame = df.filter(pl.col("group") == g)
            color: str = init.select("color").unique().collect().item()
            g_format = g.replace(" ", "_")
            out_file: Path = out_dir / f"{g_format}_heatmap.svg"
            g_df = (
                init.collect()
                .to_pandas()
                .pivot(
                    index=["chrom", "strand", "start", "group", "color"],
                    columns=self._pos_col,
                    values=self._count_col,
                )
            )
            self._plot_heatmap(
                df=g_df,
                width=width,
                height=height,
                vmin=vmin,  # type: ignore
                vmax=vmax,  # type: ignore
                cmap=self._get_cmap(color, name=g),
                group=g,
                roi_length=roi_length,
                xlabel=xlabel,
                ylabel=ylabel,
                title="Crosslink heatmap",
                outfile=out_file,
            )

    def _plot_heatmap(
        self,
        df: pd.DataFrame,
        width: float,
        height: float,
        vmin: float,
        vmax: float,
        cmap: mcolors.LinearSegmentedColormap,
        group: str,
        roi_length: Optional[int] = None,
        xlabel: str = "Relative position",
        ylabel: str = "Crosslink counts",
        title: str = "Crosslink heatmap",
        outfile: Optional[Path] = None,
    ):
        fig, ax = plt.subplots(figsize=(width, height))
        xticks = self._get_xticks(df.columns, roi_length=roi_length)
        hm = sns.heatmap(
            df,
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            yticklabels=False,
            zorder=2,
        )
        ax.set_xticks(xticks[0], labels=xticks[1])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title} for {group} samples")
        if roi_length is not None:
            roi_indices = df.columns.get_indexer_for([1, roi_length + 1])
            if not any(roi_indices < 0):
                ax.axvline(
                    roi_indices[0].item(),
                    color="black",
                    linestyle="--",
                    linewidth=1,
                    # zorder=1,
                )
                ax.axvline(
                    roi_indices[1].item(),
                    color="black",
                    linestyle="--",
                    linewidth=1,
                    # zorder=1,
                )
        for spine in ax.spines.values():
            spine.set_visible(True)
            # spine.set_linewidth(2)
            # spine.set_edgecolor('black')
        if outfile.exists():  # type: ignore
            logger.warning(
                f"Warning: {outfile} already exists and will be overwritten."
            )
        fig.savefig(str(outfile), bbox_inches="tight")
        plt.close(fig)

    def _get_xticks(
        self, columns: pd.Index, step: int = 50, roi_length: Optional[int] = None
    ) -> tuple[list[int], list[int]]:
        exact_mask = columns % step == 0
        filtered_idx = (columns[exact_mask]).to_list()
        indices = ((np.where(exact_mask)[0]) + 0.5).tolist()
        if roi_length is None:
            return indices, filtered_idx
        else:
            roi_ticks: list[int] = self._roi_ticks(roi_length)
            roi_mask = columns.isin(roi_ticks)
            if any(roi_mask):
                roi_idx = (columns[roi_mask]).to_list()
                roi_indices = ((np.where(roi_mask)[0]) + 0.5).tolist()
                return sorted(indices + roi_indices), sorted(filtered_idx + roi_idx)
            else:
                return indices, filtered_idx

    def _roi_ticks(self, roi_length: int) -> list[int]:
        """_roi_ticks Helper function
        Generate x-axis ticks for the ROI plot based on the ROI length. The ticks are generated at regular intervals of 10 bases, with the ROI region (0) in the center.
        Args:
            roi_length: (int) Length of the ROI region.

        Returns:
            list[int]: List of x-axis tick positions.
        """
        if roi_length is None:
            raise ValueError("ROI length must be provided to generate x-axis ticks.")
        if roi_length < 10:
            return [roi_length]
        if roi_length < 30:
            roi_range = roi_length // 2 * 2
            return np.arange(0, roi_range, 2).tolist()
        else:
            roi_range = roi_length // 10 * 10  # round down to nearest 10
            return np.arange(0, roi_range, 10).tolist()
