import logging
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import polars as pl
import seaborn as sns
from matplotlib import ticker
from matplotlib.lines import Line2D
from matplotlib.text import Text

logger = logging.getLogger(__name__)


class BoundedMaxNLocator(ticker.Locator):
    def __init__(self, vmin: int, vmax: int, n_ticks: int = 6):
        super().__init__()
        self.vmin = vmin
        self.vmax = vmax
        self.min_n_ticks: int = 2 if vmax > 10 else 1
        # MaxNLocator automatically handles 'nice' steps for any vmax
        self.base_locator = ticker.MaxNLocator(
            nbins=n_ticks,
            steps=[2, 5, 10],
            prune="both",
            min_n_ticks=self.min_n_ticks,
        )

    def __call__(self) -> list[int]:
        # Temporarily tell the locator the bounds to generate ticks for
        ticks = self.base_locator.tick_values(self.vmin, self.vmax)
        return [t for t in ticks if self.vmin <= t <= self.vmax]  # type: ignore


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
        ymax: float | None = None,
        errorbar: str | None = None,
        show_group_mean: bool = False,
        roi_length: int | None = None,
        labelsize: int = 9,
        xlabel: str = "Relative position",
        ylabel: str = "Counts",
        title: str = "Profile",
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
        # padding for labels
        pad = 6
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
        labels: list[int] = self._format_ticks(ax.get_xticklabels(), xmin, xmax)
        ax.set_xticks(labels)
        if roi_length is not None:
            if roi_length > 10:
                ax.xaxis.set_minor_locator(
                    BoundedMaxNLocator(vmin=0, vmax=roi_length, n_ticks=4)
                )
                minor_labelsize: int = labelsize - 1
                rotation: int = 45
                minor_pad = pad - 2
            else:
                ax.xaxis.set_minor_locator(ticker.FixedLocator([roi_length]))
                minor_labelsize = labelsize
                rotation: int = 0
                minor_pad = pad
            ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter("%d"))
            ax.tick_params(
                which="minor",
                labelsize=minor_labelsize,
                rotation=rotation,
                pad=minor_pad,
            )
            # custom ticks for ROI region if roi_length is provided
            # custom ticks for ROI region
            # vertical lines to indicate ROI region
            ax.axvline(0, color="black", linestyle="--", linewidth=1, zorder=0)
            ax.axvline(
                roi_length + 1, color="black", linestyle="--", linewidth=1, zorder=0
            )

        # set y limit if provided
        if ymax is not None:
            ax.set_ylim(top=ymax)
        # Labels
        ax.tick_params(
            which="major",
            labelsize=labelsize,
            pad=pad,
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        # legend
        lines, clabels = self._per_group_legend(self._group_cmap)
        ax.legend(lines, clabels, title=None, loc="upper right", bbox_to_anchor=(1, 1))
        # save
        if Path(output).exists():
            logger.warning(f"Warning: {output} already exists and will be overwritten.")
        fig.savefig(output, bbox_inches="tight")
        plt.close(fig)

    def _format_ticks(self, ticks: list[Text], xmin: int, xmax: int) -> list[int]:
        """_format_ticks Helper function

        Format x-axis tick labels as integers and filter them based on the specified range.

        Args:
            ticks: (list[Text]) List of x-axis tick labels as Text objects.
            xmin: (int) Minimum x-axis value to include.
            xmax: (int) Maximum x-axis value to include.

        Returns:
            list[int]: List of formatted and filtered x-axis tick positions.
        """
        labels: list[int] = []
        for tick in ticks:
            tlabel = re.sub(
                r"[\u2212\u2010\u2011\u2012\u2013\u2014\u2015]",
                "-",
                tick.get_text(),
            )
            ilabel = int(tlabel)
            if (ilabel < xmin) or (ilabel > xmax):
                continue
            labels.append(ilabel)
        return labels

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
