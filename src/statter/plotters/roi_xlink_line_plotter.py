import numpy as np
import pandas as pd
import seaborn as sns


class LinePlotter:

    def __init__(
        self, df: pd.DataFrame, meta_df: pd.DataFrame, smoothing_window: int = 5
    ) -> None:
        self.df = self._grouper_smoother(df, smoothing_window)
        self.meta_df = meta_df

    def _grouper_smoother(
        self, df: pd.DataFrame, smoothing_window: int
    ) -> pd.DataFrame:
        """_grouper_smoother Helper function
        Group dataframe by group and sample, then apply smoothing to the values. The smoothing is done by convolving the values with a kernel.
        Args:
            df: (pd.DataFrame) DataFrame to be grouped and smoothed. Must contain 'group' and 'sample' columns.
            smoothing_window: (int) Size of the smoothing window. If 2, a simple average of the current and next value is taken. If greater than 2, a Gaussian kernel is used for smoothing.

        Returns:
            pd.DataFrame: The grouped and smoothed DataFrame.
        """

        def _smoother(row, kernel):
            return np.convolve(row, kernel, mode="same")

        req_cols: set[str] = {"group", "sample"}
        if not req_cols.issubset(set(df.columns)):
            raise ValueError(
                f"DataFrame must contain the following columns: {', '.join(req_cols)}"
            )
        df = df.groupby(list(req_cols)).mean()
        if smoothing_window < 2:
            return df
        kernel = self._kernel(smoothing_window)
        old_cols = df.columns
        df = df.apply(
            lambda row: _smoother(row, kernel), axis=1, result_type="expand"
        )  # type: ignore
        df.columns = old_cols
        return df

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

    def plot(self):
        var_cols: pd.Index = self.df.columns
        id_cols: pd.api.typing.FrozenList = self.df.index.names
        df_melt: pd.DataFrame = self.df.reset_index().melt(
            id_vars=id_cols,
            value_vars=var_cols,
            var_name="position",
            value_name="count",
        )
        lp = sns.lineplot(
            data=df_melt,
            x="position",
            y="count",
            hue="group",
            units="sample",
            estimator=None,
        )


"""
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Create the figure and axes objects
fig, ax = plt.subplots(figsize=(10, 6))

# 2. Draw the individual sample lines (Background layer)
sns.lineplot(
    data=df_agg, x='position', y='count', hue='group', 
    units='sample', estimator=None, 
    alpha=0.3, legend=False, ax=ax
)

# 3. Draw the overall group mean (Foreground layer)
sns.lineplot(
    data=df_agg, x='position', y='count', hue='group', 
    estimator='mean', linewidth=4, ax=ax
)

# 4. Set the Y-axis limits
# You can provide (min, max). Use None for one side to let it auto-scale.
ax.set_ylim(0, 500) 

# Optional: Add a title or clean up labels
ax.set_title('Count by Position per Group (Individual Samples vs. Mean)')
plt.show()

add custom legends

from matplotlib.lines import Line2D

fig, ax = plt.subplots(figsize=(10, 6))

# Define a palette so colors are consistent
my_palette = {'A': 'royalblue', 'B': 'firebrick'}

sns.lineplot(
    data=df_agg, x='position', y='count', 
    hue='group', units='sample', estimator=None, 
    palette=my_palette, alpha=0.6, ax=ax, legend=False # Turn off auto-legend
)

# 1. Create custom "handles" (the icons in the legend)
custom_lines = [
    Line2D([0], [0], color='royalblue', lw=2),
    Line2D([0], [0], color='firebrick', lw=2)
]

# 2. Apply the custom legend
ax.legend(custom_lines, ['Group A (Normal)', 'Group B (Variant)'], 
          title="Legend", loc='upper right', frameon=False)

ax.set_ylim(0, 500)
plt.show()

"""
