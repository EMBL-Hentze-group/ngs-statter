import gzip
import logging
from functools import partial
from io import BufferedReader
from itertools import islice
from pathlib import Path
from typing import Callable
import json

"""
parse fastq file and write read length stats to a tsv file
"""


class FqLength:
    """
    Parse a given fastq file and write length stats to a tsv file
    """

    def __init__(self, fq: str, out_file: str) -> None:
        """
        fq: /path/to/fastq.file, either plain text or gzipped
        out_file: /path/to/stats.csv file
        sample_name: Optional, sample name for this fastq file. If not given, base file name will be used as sample name
        """
        self.fq = fq
        self.out_file = out_file
        self._fh = self._file_reader()
        self._read_len_stats = {}

    def _file_reader(self) -> Callable:
        """
        Helper function, return appropriate file handler, io stream and decoder
        """
        with open(self.fq, "rb") as _rb:
            if _rb.read(2) == b"\x1f\x8b":
                return partial(gzip.open, mode="rb")
            else:
                return partial(open, mode="rb")

    def read_length(self) -> None:
        """
        Compute and write read length stats
        """
        with BufferedReader(self._fh(self.fq)) as fh:
            while True:
                seq = list(islice(fh, 4))
                if not seq:
                    break
                elif len(seq) < 4 or seq[0].decode("utf-8")[0] != "@":
                    continue
                read_length = len(seq[1].decode("utf-8").strip())
                try:
                    self._read_len_stats[read_length] += 1
                except KeyError:
                    self._read_len_stats[read_length] = 1
        self._write_stats()

    def _write_stats(self) -> None:
        """
        Helper function
        Write read length stats to a file
        """
        if Path(self.out_file).exists():
            logging.warning(f"Re-writing file {self.out_file}")
        with open(self.out_file, "w") as wh:
            # sort and dump read length dictionary
            json.dump(dict(sorted(self._read_len_stats.items())), wh)
