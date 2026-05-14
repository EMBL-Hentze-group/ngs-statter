import logging
import re
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
from matplotlib.lines import Line2D

logger = logging.getLogger(__name__)


class LinePlotter:

    def __init__(self, df_path: Path | str) -> None:
        self._req_cols: list[str] = [
            "chrom",
            "start",
            "strand",
            "group",
            "sample",
            "color",
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

    @property
    def _group_cmap(self) -> dict[str, str]:
        df_cols: pl.DataFrame = self._df.select(["group", "color"]).unique().collect()
        return dict(zip(df_cols["group"], df_cols["color"]))

    @property
    def _sample_cmap(self) -> dict[str, str]:
        df_cols: pl.DataFrame = self._df.select(["sample", "color"]).unique().collect()
        return dict(zip(df_cols["sample"], df_cols["color"]))

    def plot(
        self,
        output: str,
        width: float = 10,
        height: float = 6,
        ymax: Optional[float] = None,
        errorbar: Optional[str] = None,
        show_group_mean: bool = False,
        roi_length: Optional[int] = None,
        xlabel: str = "Relative position",
        ylabel: str = "Crosslink counts",
        title: str = "Crosslink profile",
    ) -> None:
        df_melt: pd.DataFrame = (
            self._df.unpivot(
                self._rel_pos,
                index=self._req_cols,
                variable_name="position",
                value_name="count",
            )
            .with_columns(pl.col("position").cast(pl.Int32))
            .collect()
            .to_pandas()
        )
        # plotting
        fig, ax = plt.subplots(figsize=(width, height))
        # per sample
        per_sample = sns.lineplot(
            data=df_melt,
            x="position",
            y="count",
            hue="sample",
            estimator="mean",
            errorbar=errorbar,
            legend=False,
            alpha=0.40 if show_group_mean else 0.95,
            palette=self._sample_cmap,
            zorder=2,
            ax=ax,
        )
        if show_group_mean:
            # per group mean
            per_group = sns.lineplot(
                data=df_melt,
                x="position",
                y="count",
                hue="group",
                estimator="mean",
                errorbar=None,
                legend=False,
                linewidth=3.0,
                alpha=0.80,
                palette=self._group_cmap,
                zorder=1,
                ax=ax,
            )
        # set x-axis limits to nearest 10
        xmin = (df_melt["position"].min() // 10 - 1) * 10
        xmax = (df_melt["position"].max() // 10 + 1) * 10
        ax.set_xlim(xmin, xmax)
        if roi_length is not None:
            # custom ticks for ROI region if roi_length is provided
            main_ticks: list[int] = np.arange(
                (xmin // 50 + 1) * 50, (xmax // 50 + 1) * 50, 50
            ).tolist()
            ticks: list[int] = sorted(main_ticks + self._roi_ticks(roi_length))
            # custom ticks for ROI region
            ax.set_xticks(ticks)
            # vertical lines to indicate ROI region
            ax.axvline(0, color="black", linestyle="--", linewidth=1, zorder=0)
            ax.axvline(
                roi_length + 1, color="black", linestyle="--", linewidth=1, zorder=0
            )
        # set y limit if provided
        if ymax is not None:
            ax.set_ylim(top=ymax)
        # Labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        # legend
        lines, labels = self._per_group_legend(self._group_cmap)
        ax.legend(lines, labels, title=None, loc="upper right", bbox_to_anchor=(1, 1))
        # save
        if Path(output).exists():
            logger.warning(f"Warning: {output} already exists and will be overwritten.")
        fig.savefig(output, bbox_inches="tight")
        plt.close(fig)

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
            return [roi_length // 2, roi_length]
        if roi_length < 30:
            roi_range = roi_length // 2 * 2
            return np.arange(0, roi_range, 2).tolist()
        else:
            roi_range = roi_length // 10 * 10  # round down to nearest 10
            return np.arange(0, roi_range, 10).tolist()

    def _per_group_legend(self, cmap: dict[str, str]) -> tuple[list[Line2D], list[str]]:
        """
        _per_group_legend Helper function
        Generate custom legend handles and labels for the per-group plot.
        Args:
            cmap: (dict[str, str]) Dictionary mapping group names to colors.

        Returns:
            tuple[list[Line2D], list[str]]: Tuple containing a list of Line2D handles and a list of group names.
        """
        items: list[str] = sorted(cmap.keys())
        colors: list[str] = [cmap[k] for k in items]
        handles: list[Line2D] = [Line2D([0], [0], color=c, lw=2) for c in colors]
        return (handles, items)
