import gzip
import json
import logging
from collections import defaultdict
from functools import partial
from io import BufferedReader
from itertools import islice
from pathlib import Path
from typing import Callable

"""
Read unclassified read ids from ganon2 output, read the corresponding fastq file and write ganon2 unclassified reads into a separate fastq file
"""
logger = logging.getLogger(__name__)


class GanonUnclassifiedFq:
    """
    Write Ganon2 unclassified reads into a separate fastq file
    """

    def __init__(
        self, unc: str, src_fq: str, target_fq: str, read_length_json: str
    ) -> None:
        """
        Arguments:
         unc: ganon unclassified read ids file
         src_fq: source fastq file
         target_fq: target fastq file
         read_length_json: target fastq read length distribution
        """
        self.unc = unc
        self.src_fq = src_fq
        self.target_fq = target_fq
        self.read_length_json = read_length_json
        self.fastq_unc_ids = set()  # unclassified fastq ids
        self._unc_reader()
        self._read_lengths = defaultdict(dict)  # read length dictionary

    def _unc_reader(self) -> None:
        """
        Read ganon2 unclassified fastq read ids
        """
        with open(self.unc, "rb") as _rb:
            if _rb.read(2) == b"\x1f\x8b":
                logger.debug(f"{self.unc} is gzipped")
                unc_reader = partial(gzip.open, mode="rt")
            else:
                logger.debug(f"{self.unc} is plain text")
                unc_reader = partial(open, mode="r")
        with unc_reader(self.unc) as _uh:
            for read_id in _uh:
                self.fastq_unc_ids.add(read_id.strip())
        if len(self.fastq_unc_ids) == 0:
            raise RuntimeError(f"Cannot parse fastq read ids from {self.unc}")
        logger.info(f"Found {len(self.fastq_unc_ids)} read ids in {self.unc} ")

    @staticmethod
    def _file_reader(file_name) -> Callable:
        """
        Helper function, return appropriate file handler, io stream and decoder
        """
        with open(file_name, "rb") as _rb:
            if _rb.read(2) == b"\x1f\x8b":
                logger.debug(f"{file_name} is gzipped")
                return partial(gzip.open, mode="rb")
            else:
                logger.debug(f"{file_name} is plain text")
                return partial(open, mode="rb")

    @staticmethod
    def _file_writer(fq) -> Callable:
        """
        Helper function
        return appropriate file handler based on suffix
        """
        fqp = Path(fq)
        if fqp.exists():
            logger.warning(f"Re-writing file {fq}")
        fq_suffix = fqp.suffixes[-1].lower()
        if fq_suffix == ".gz":
            return partial(gzip.open, mode="wb", compresslevel=6)
        else:
            return partial(open, mode="wb")

    def write_unclassified_reads(self):
        """
        Write ganon2 unclassified reads into a separate fastq file
        """
        fq_reader = self._file_reader(self.src_fq)
        fq_writer = self._file_writer(self.target_fq)
        with BufferedReader(fq_reader(self.src_fq)) as freader, fq_writer(
            self.target_fq
        ) as fwriter:
            while True:
                seq = list(islice(freader, 4))
                if not seq:
                    break
                elif len(seq) < 4:
                    continue
                header = seq[0].decode("utf-8").strip()
                read_length = len(seq[1].decode("utf-8").strip())
                header_id = header[1:].split(" ")
                if header[0] != "@":
                    continue
                elif header_id[0] not in self.fastq_unc_ids:
                    self._add_read_length_count("contamination", read_length)
                    continue
                self._add_read_length_count("unclassified", read_length)
                fwriter.write(b"".join(seq))
        self._write_read_length_json()

    def _add_read_length_count(self, read_type: str, read_length: int) -> None:
        """
        Helper function
        add read length count
        Arguments:
            read_type: string, either "contamination" or "unclassified"
            read_length: int, read length
        """
        try:
            self._read_lengths[read_type][read_length] += 1
        except KeyError:
            self._read_lengths[read_type][read_length] = 1

    def _write_read_length_json(self):
        """
        Helper function
        Write read length distribution in json for ganon2 unclassified reads
        """
        if Path(self.read_length_json).exists():
            logger.warning(f"Re-writing {self.read_length_json}")
        with open(self.read_length_json, "w") as jh:
            json.dump(self._read_lengths, jh)
