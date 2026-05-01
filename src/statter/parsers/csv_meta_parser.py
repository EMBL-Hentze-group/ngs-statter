import logging
import sys
from pathlib import Path

import pandas as pd
from matplotlib import colors

logger = logging.getLogger(__name__)


class MetaReader:
    """
    Check sanity of the supplied yaml file
    """

    def __init__(self, metafile: str | Path) -> None:
        # each file must have the following keys
        # name: name to call the sample
        # group: group to which the sample belongs
        self._req_cols: set[str] = {"file", "sample", "group"}
        self.opt_cols: set[str] = {"color"}
        self.metafile = metafile

    def read_meta(self) -> pd.DataFrame:
        """meta_reader read the metadata file and check if it has the required columns

        Raises:
            ValueError: if the metadata file does not have the required columns
        Returns:
            pd.DataFrame: DataFrame containing the metadata information
        """
        metadf: pd.DataFrame = pd.read_csv(self.metafile, sep="\t")
        missing_cols = self._req_cols - set(metadf.columns)
        if missing_cols:
            raise ValueError(
                f"Metadata file {self.metafile} is missing required columns: {', '.join(missing_cols)}"
            )
        return metadf

    def yaml_example(self):
        example = """
        A compatible yaml file should like the following:
        
        $ cat example.yaml
        /path/to/input_1_file.bed(.gz): # ":" MUST be there, shoji/htseq-clip "sites" file per sample
            name: IP1 #  a suitable name for the sample
            group: IP # group to which sample belongs
        /path/to/input_2_file.bed: # 
            name: SMI1 # a suitable name for the sample
            group: SMI
        /path/to/input_3_file.bed: # 
            name: IP2 # a suitable name for the sample
            group: IP
        
        FYI: "name:" and "group:" lines must begin with a single space character.
        """
        sys.stdout.write(example + "\n")

    def metadata_example(self) -> None:
        example = """
        A compatible metadata file should like the following:
        
        $ cat example_metadata.csv
        file<\t>sample<\t>group<\t>color[optional]
        /path/to/input_1_file.bed(.gz)<\t>IP1<\t>IP<\t>blue
        /path/to/input_2_file.bed<\t>SMI1<\t>SMI<\t>#FF0000
        /path/to/input_3_file.bed<\t>IP2<\t>IP<\t>green

        FYI: columns should be separated by <tab> character
            the first line (header) should contain the column names "file", "sample", "group" and optionally "color"
            color can either be a color name (e.g. "blue", "red", "green") or a hex code (e.g. "#FF0000")
        """
        sys.stdout.write(example + "\n")
