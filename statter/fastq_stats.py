import gzip
import logging
import re
import gzip
from functools import partial
from io import BufferedReader
from itertools import islice
from pathlib import Path
from typing import Callable, Optional

from tqdm import tqdm

"""
parse fastq file and write read length stats to a tsv file
"""


class FqLength:
    """
    Parse a given fastq file and write length stats to a tsv file
    """

    def __init__(
        self, fq: str, out_file: str, sample_name: Optional[str] = None
    ) -> None:
        """
        fq: /path/to/fastq.file, either plain text or gzipped
        out_file: /path/to/stats.csv file
        sample_name: Optional, sample name for this fastq file. If not given, base file name will be used as sample name
        """
        self.fq = fq
        self.out_file = out_file
        if sample_name is None:
            self.sample_name = re.sub("\..*$", "", Path(self.fq).stem)
        else:
            self.sample_name = sample_name
        logging.info(f"Using {self.sample_name} as sample name")
