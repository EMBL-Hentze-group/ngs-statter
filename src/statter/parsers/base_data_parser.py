import logging
import tempfile
from abc import ABC, abstractmethod
from collections.abc import Callable
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)


class BaseDataFrameProcessor(ABC):
    """Shared pipeline for metadata-driven tabular genomic parsers."""

    def __init__(
        self,
        metadf: pl.DataFrame,
        tempfolder: str | Path | None,
    ) -> None:
        self.metadata = metadf

        self._tmpdir: Path = self._init_tmpdir(tempfolder)

        self._groups: list[str] = (
            self.metadata.get_column("group").unique().sort().to_list()
        )
        self._samples: list[str] = (
            self.metadata.get_column("sample").unique().sort().to_list()
        )

    def _init_tmpdir(self, tempfolder: str | Path | None) -> Path:
        """_init_tmpdir Helper function
        Generate a temporary directory for intermediate file storage.
        Args:
            tempfolder: Optional base folder for the temporary directory. If None, the system temporary directory is used.

        Returns:
            Path to the created temporary directory.
        """
        if tempfolder is None:
            return Path(tempfile.gettempdir()) / next(tempfile._get_candidate_names())  # type: ignore
        tmpdir = Path(tempfolder) / next(tempfile._get_candidate_names())  # type: ignore
        tmpdir.mkdir(exist_ok=True, parents=True)
        return tmpdir

    def _get_tmp_fname(self, suffix: str | None = None) -> Path:
        fname = next(tempfile._get_candidate_names())  # type: ignore
        if suffix is not None:
            fname = f"{fname}{suffix}"
        return self._tmpdir / fname

    @property
    def tmpdir(self) -> Path:
        return self._tmpdir

    @property
    @abstractmethod
    def _schema(self) -> dict[str, pl.DataType]:
        """CSV/TSV schema for the input file."""

    @property
    @abstractmethod
    def _partition_cols(self) -> list[str]:
        """Columns used to partition parquet output."""

    @property
    def partition_cols(self) -> list[str]:
        return self._partition_cols

    @property
    @abstractmethod
    def _group_by_cols(self) -> list[str]:
        """Columns used to group data before processing."""

    @property
    @abstractmethod
    def _rename_cols(self) -> dict[str, str] | None:
        """Columns to rename in the input data."""

    @abstractmethod
    def _parse_and_transform(
        self, filepath: str, sample: str, group: str
    ) -> pl.LazyFrame:
        """Parse one file and return a transformed lazy frame ready for write-out."""

    @staticmethod
    def _temp_name_provider(args) -> str:
        """Generate partitioned parquet output filenames."""
        fps: list[str] = []
        df = args.partition_keys
        for col in df.columns:
            fps.append(f"{col}={df[col][0]}")
        path = Path().joinpath(*fps) / f"{next(tempfile._get_candidate_names())}.parquet"  # type: ignore
        return str(path)

    def to_signal(self) -> Path:
        """Iterate metadata rows, parse each file, and write partitioned parquet output."""
        out_folder: Path = self._tmpdir / next(tempfile._get_candidate_names())  # type: ignore
        out_folder.mkdir(exist_ok=True, parents=True)

        for dat in self.metadata.iter_rows(named=True):
            logger.info(
                f"Processing file {dat['file']} for sample {dat['sample']} in group {dat['group']}"
            )
            lf = self._parse_and_transform(
                filepath=dat["file"],
                sample=dat["sample"],
                group=dat["group"],
            )
            lf.sink_parquet(
                pl.PartitionBy(
                    base_path=out_folder,
                    file_path_provider=self._temp_name_provider,
                    key=self._partition_cols,
                    include_key=False,
                ),
                mkdir=True,
            )

        return out_folder
