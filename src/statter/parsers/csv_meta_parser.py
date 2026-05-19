import logging
from csv import DictReader
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from matplotlib import colors
from numpy import linspace
from rich.console import Console
from rich.table import Table

from statter.json_models.sample_meta import SampleMeta

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
        with open(self.metafile, "r") as fh:
            reader = DictReader(fh, delimiter="\t")
            missing = SampleMeta.required_fields() - set(reader.fieldnames)  # type: ignore
            if missing:
                raise ValueError(
                    f"Metadata file {self.metafile} is missing required columns: '{', '.join(missing)}'. Please check the metadata file!"
                )
            self._metadf = pl.DataFrame(
                [SampleMeta.model_validate(row) for row in reader]
            )

        self._check_groups()
        self._check_samples()
        self._check_colors()
        return self._metadf

    def _check_groups(self) -> None:
        """
        Check that samples >= groups
        """
        if self._metadf["sample"].unique().len() < self._metadf["group"].unique().len():
            raise ValueError(
                f"Metadata file {self.metafile} contains more unique groups than samples. Please check the metadata file and ensure that each sample is assigned to a valid group."
            )

    def _check_samples(self) -> None:
        """
        Check that samples are unique
        """
        if self._metadf["sample"].unique().len() != self._metadf.height:
            raise ValueError(
                f"Metadata file {self.metafile} contains duplicate sample names. Please check the metadata file and ensure that each sample has a unique name."
            )

    def _check_colors(self) -> None:
        """_check_colors Helper function
        Checks if the metadata DataFrame contains a valid "color" column. If not, generates colors randomly
        """
        if any(self._metadf["color"].is_null()):
            self._generate_colors()
        elif (
            self._metadf["group"].unique().len() != self._metadf["color"].unique().len()
        ):
            logger.warning(
                f"Metadata file {self.metafile} contains either invalid color values or mismatched number of colors and groups. Colors will be generated randomly."
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
        console = Console()
        table = Table(title="Example Metadata File")

        table.add_column("file")
        table.add_column("sample")
        table.add_column("group")
        table.add_column(
            "color *optional*",
        )

        table.add_row("/path/to/input_1_file.bed(.gz)", "IP1", "IP", "blue")
        table.add_row("/path/to/input_2_file.bed", "SMI1", "SMI", "#FF0000")
        table.add_row("/path/to/input_3_file.bed", "IP2", "IP", "blue")

        console.print(table)
        console.print(
            "[bold]FYI:[/bold] \n"
            "- columns should be separated by [bold]<tab>[/bold](\\t) character\n"
            "- The first line (header) should contain the column names 'file', 'sample', 'group' and optionally 'color'.\n"
            "- The order of columns does not matter, but the column names are case sensitive.\n\n"
            "Color can either be a color name (e.g. 'blue', 'red', 'green') or a hex code (e.g. '#FF0000').\n"
            "Colors are per group, so if multiple samples belong to the same group, they should have the same color. "
            "If colors are not provided or invalid, they will be generated randomly."
        )


def main():
    # meta_path1 = "test_data/crosslinks/soniCLIP_Histone_3UTR_crosslinks_naln_5.csv"
    # meta_reader1 = MetaReader(meta_path1)
    # print(meta_reader1.read_meta())

    # meta_path2 = "test_data/crosslinks/jumbled.csv"
    # meta_reader2 = MetaReader(meta_path2)
    # print(meta_reader2.read_meta())
    print(MetaReader.metadata_example())


if __name__ == "__main__":
    main()
