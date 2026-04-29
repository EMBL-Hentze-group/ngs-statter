import gzip
import logging
from pathlib import Path
from typing import List, TextIO

from pyfastx import Fastx

logger = logging.getLogger(__name__)


class HeaderFixer:
    """
    Fix UMI positions in FASTQ headers produced by Flexbar.

    Flexbar appends UMIs to the end of the headers, which means that in Fastq files from Casava >1.8 versions
    the UMI is appended to the end of the comment section. All aligners typically remove this part of the header
    during alignment, which results in loss of UMI information.

    This class moves the UMI to the enf of the actual read identifier, before the first space separator.
    """

    def __init__(self, in_fq: str, out_fq: str, separator: str = "_") -> None:
        self.in_fq = in_fq
        self.out_fq = out_fq
        self.separator = separator

    def _writer(self) -> TextIO:
        out: Path = Path(self.out_fq)
        if out.suffix in set([".gz", ".gzip"]):
            return gzip.open(self.out_fq, "wt")
        else:
            return open(self.out_fq, "w")

    def fix_umi_header(self) -> None:
        """
        Process the input FASTQ file and write the modified headers to the output FASTQ file.
        """
        warnings: List[str] = []
        n_warnigs: int = 0
        buffer: List[str] = []
        buf_size: int = 50000
        with self._writer() as _fh:
            for name, seq, qual, comment in Fastx(self.in_fq, comment=True):
                if comment is None:
                    if len(warnings) < 10:
                        warnings.append("Detected header without comment section.")
                        n_warnigs += 1
                    buffer.append(f"@{name}\n{seq}\n+\n{qual}\n")
                    continue
                dat: List[str] = comment.split(self.separator)
                if len(dat) < 2:
                    if len(warnings) < 10:
                        warnings.append(
                            f"Header '{name}' does not contain expected UMI separator '{self.separator}'"
                        )
                        n_warnigs += 1
                    buffer.append(f"@{name} {comment}\n{seq}\n+\n{qual}\n")
                    continue
                umi: str = dat[-1]
                new_header: str = f"{name}{self.separator}{umi}"
                buffer.append(f"@{new_header} {' '.join(dat[:-1])}\n{seq}\n+\n{qual}\n")
                if len(buffer) >= buf_size:
                    _fh.writelines(buffer)
                    buffer.clear()
            if buffer:
                _fh.writelines(buffer)
        if warnings:
            for w in warnings:
                logger.warning(w)
            if n_warnigs > len(warnings):
                logger.warning(f"... and {n_warnigs - len(warnings)} more warnings.")
