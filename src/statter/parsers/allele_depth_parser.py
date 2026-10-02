import logging
from pathlib import Path
from warnings import warn

import polars as pl
from pysam import bcftools

from statter.parsers.base_data_parser import BaseDataFrameProcessor

logger = logging.getLogger(__name__)


class AlleleDepthParser(BaseDataFrameProcessor):
    """AlleleDepthParser
    Parse Allele depth data from Variant calling and return minor allele read depth.
    The parser expects a CSV file with columns: CHROM, POS, REF, ALT, sample, AD0, AD1. from biallelic vcf files.
    If there are multiple alternate alleles, read depth for all alternate alleles per CHROM and POS combinations are summed up

    Args:
        BaseDataFrameProcessor: _description_
    """

    def __init__(
        self,
        metadf: pl.DataFrame,
        tempfolder: str | Path | None,
        genome_fa: str | Path | None,
        variant_type: str | None,
        dp: int = 0,
        ncores: int = 1,
    ):
        super().__init__(metadf, tempfolder=tempfolder)
        self.genome_fa = genome_fa
        self.variant_type = variant_type
        self.dp = dp
        self.ncores = ncores
        self._variant_filter: dict[str, str] = {
            "deletion": "'indel' && strlen(REF) > strlen(ALT)",
            "insertion": "'indel' && strlen(REF) < strlen(ALT)",
        }

    def _filter_builder(self) -> str:
        """_filter_builder Helper function
        Construct filter strings for variant types and depth thresholds. The filter string can be used in bcftools query to filter variants based on type and depth.

        Returns:
            filter string for bcftools query based on the specified variant type and depth threshold.
        """
        filters = []
        if self.variant_type is not None:
            variant_type = self.variant_type.lower()
            if variant_type in self._variant_filter:
                filters.append(f"TYPE={self._variant_filter[variant_type]}")
            else:
                filters.append(f"TYPE='{variant_type}'")
        if self.dp > 0:
            filters.append(f"FMT/DP>={self.dp}")
        return " && ".join(filters)

    @property
    def _schema(self) -> dict[str, pl.DataType]:
        return {
            "chrom": pl.String,
            "end": pl.Int64,
            "sample": pl.String,
            "count": pl.String,
        }  # type: ignore

    @property
    def _group_by_cols(self) -> list[str]:
        return ["chrom", "end"]

    @property
    def _partition_cols(self) -> list[str]:
        # intentionally in small letters to match the renamed columns after parsing
        return ["chrom"]

    @property
    def _rename_cols(self) -> dict[str, str] | None:
        return None

    def _norm(self, filepath: str, sample: str):
        """_norm Helper function
        Call bcftools norm to normalize variants, if genome is given,
        convert all multiallelic variants to biallelic for easy processing
        Args:
            filepath: path to the input VCF file.
            sample: sample name corresponding to the VCF file.
        Returns:
            path to the normalized VCF file.
        """
        _tmp_vcf: str = str(self._get_tmp_fname(suffix=f"{sample}_.vcf.gz"))
        norm_opts: list[str] = [
            "--threads",
            str(self.ncores),
            "-m",
            "-any",
            "-Wtbi",
            "-Oz",
            "-o",
            _tmp_vcf,
        ]
        if self.genome_fa:
            norm_opts.extend(["-f", str(self.genome_fa)])
        norm_opts.append(filepath)  # input VCF at the end
        bcftools.norm(*norm_opts, catch_stdout=False)
        return _tmp_vcf

    def _query(self, vcf: str, sample: str):
        """_query Helper function
        Call bcftools query to extract allele depth information from the VCF file.
        Args:
            vcf: path to the input VCF file.
            sample: sample name corresponding to the VCF file.

        Returns:
            path to the CSV file containing the extracted allele depth information.
        """
        _tmp_csv = str(self._get_tmp_fname(suffix=f"{sample}_.csv"))
        query_opts: list[str] = [
            "-f",
            "%CHROM\t%POS\t[%SAMPLE\t%AD{1}\n]",
            "-o",
            _tmp_csv,
        ]
        query_filter = self._filter_builder()
        if len(query_filter) > 0:
            query_opts.extend(["-i", query_filter])
        query_opts.append(vcf)  # input VCF at the end
        bcftools.query(*query_opts, catch_stdout=False)
        return _tmp_csv

    def _parse_and_transform(
        self, filepath: str, sample: str, group: str
    ) -> pl.LazyFrame:
        # bcftools norm and query to preprocess the VCF file
        logger.info(f"Processing file {filepath} : sample {sample},  group {group}")
        query_csv: str = self._query(self._norm(filepath, sample), sample)
        # process output
        ad = pl.scan_csv(
            query_csv,
            has_header=False,
            separator="\t",
            schema=self._schema,
        )
        ad_sample: list[str] = (
            ad.select(pl.col("sample"))
            .unique()
            .collect()
            .get_column("sample")
            .to_list()
        )
        if len(ad_sample) != 1:
            raise ValueError(
                f"Expected exactly one unique sample in the AD file, but found {len(ad_sample)}: {ad_sample}"
            )
        if ad_sample[0] != sample:
            warn(
                f"Sample in the AD file ({ad_sample[0]}) does not match the sample name from metadata ({sample})!",
                RuntimeWarning,
            )
        ad = (
            ad.select(["chrom", "end", "count"])
            .with_columns(
                (pl.col("count").str.replace(".", "0").cast(pl.Int64)).alias("count")
            )
            .group_by(self._group_by_cols)
            .agg(pl.sum("count").alias("count"))
            .with_columns(
                pl.lit(sample).alias("sample"),
                pl.lit(group).alias("group"),
            )
        )
        return ad.with_columns(pl.lit("*").alias("strand"))
