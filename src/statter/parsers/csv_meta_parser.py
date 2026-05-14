import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from matplotlib import colors
from numpy import linspace

logger = logging.getLogger(__name__)


class MetaReader:
    """
    Reads and validates a metadata TSV file using polars, enforcing required columns and schema.
    """

    def __init__(self, metafile: str | Path) -> None:
        """__init__
        Read and validate the metadata file, ensuring required columns are present and colors are valid or generated.

        Args:
            metafile: (str | Path) Path to the metadata file
        """
        self._req_cols: set[str] = {"file", "sample", "group"}
        self._opt_cols: set[str] = {"color"}
        self.metafile = metafile
        # Define the schema for polars (all as string, color is optional)
        self._metadf: pl.DataFrame = pl.DataFrame()

    def read_meta(self) -> pl.DataFrame:
        """
        Reads the metadata file as a polars DataFrame and checks for required columns.

        Raises:
            ValueError: if the metadata file does not have the required columns
        Returns:
            pl.DataFrame: DataFrame containing the metadata information
        """
        # Read with relaxed schema (color is optional)
        self._metadf = pl.read_csv(
            self.metafile,
            separator="\t",
            has_header=True,
            # schema=self.schema,
            ignore_errors=True,
        )
        missing_cols = set(self._req_cols) - set(self._metadf.columns)
        if missing_cols:
            raise ValueError(
                f"Metadata file {self.metafile} is missing required columns: {', '.join(missing_cols)}"
            )
        self._check_colors()
        # check how
        return self._metadf

    def _check_colors(self) -> None:
        """_check_colors Helper function
        Checks if the metadata DataFrame contains a valid "color" column. If not, generates colors randomly
        """
        if "color" not in self._metadf.columns:
            logger.warning(
                f"Metadata file {self.metafile} is missing colors for either some or all rows. Colors will be generated automatically."
            )
            self._generate_colors()
        else:
            # validate colors
            color_check: bool = (
                self._metadf["color"]
                .map_elements(
                    lambda c: colors.is_color_like(c), return_dtype=pl.Boolean
                )
                .all()
            )
            grp_color = (
                self._metadf["group"].unique().len()
                == self._metadf["color"].unique().len()
            )  # check if number of unique colors matches number of unique groups
            if (not color_check) or (not grp_color):
                logger.warning(
                    f"Metadata file {self.metafile} contains either invalid color values or mismatched number of colors and groups. Colors will be generated automatically."
                )
                self._generate_colors()

    def _generate_colors(self) -> None:
        """_generate_colors Helper function
        Generates colors for each unique group in the metadata DataFrame using a colormap.
        Args:
            metadf: (pl.DataFrame) DataFrame containing the metadata information

        Returns:
            pl.DataFrame: DataFrame with generated colors
        """
        self._metadf = self._metadf.drop("color", strict=False)
        ncol: int = self._metadf["group"].unique().len()
        col_df: pl.DataFrame = pl.DataFrame(
            {
                "group": self._metadf["group"].unique(),
                "color": [
                    colors.to_hex(icol)
                    for icol in plt.colormaps["viridis"](linspace(0, 1, ncol))
                ],
            }
        )
        self._metadf = self._metadf.join(col_df, on="group", how="left")

    @property
    def per_group_colors(self) -> dict[str, str]:
        """per_group_colors Helper function
        Get a dictionary mapping each group to its assigned color.

        Returns:
            dict[str, str]: A dictionary where keys are group names and values are color codes.
        """
        if self._metadf is None:
            raise ValueError("Metadata has not been read yet. Call read_meta() first.")
        return dict(zip(self._metadf["group"], self._metadf["color"]))

    @property
    def per_sample_colors(self) -> dict[str, str]:
        """per_sample_colors Helper function
        Get a dictionary mapping each sample to its assigned color based on its group.

        Returns:
            dict[str, str]: A dictionary where keys are sample names and values are color codes.
        """
        if self._metadf is None:
            raise ValueError("Metadata has not been read yet. Call read_meta() first.")
        return dict(zip(self._metadf["sample"], self._metadf["color"]))

    @staticmethod
    def metadata_example() -> None:
        example = """
        A compatible metadata file should look like the following:
        
        $ cat example_metadata.csv
        file\tsample\tgroup\tcolor[optional]
        /path/to/input_1_file.bed(.gz)\tIP1\tIP\tblue
        /path/to/input_2_file.bed\tSMI1\tSMI\t#FF0000
        /path/to/input_3_file.bed\tIP2\tIP\tblue

        FYI: columns should be separated by <tab> character

        the first line (header) should contain the column names "file", "sample", "group" and optionally "color"
        color can either be a color name (e.g. "blue", "red", "green") or a hex code (e.g. "#FF0000")

        Colors are per group, so if multiple samples belong to the same group, they should have the same color. If colors are not provided or invalid, they will be generated automatically.
        """
        sys.stdout.write(example + "\n")
