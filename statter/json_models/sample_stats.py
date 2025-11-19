from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


class BaseStats(BaseModel):
    """BaseStats
    Pydantic model for basic sample statistics
    """

    sample_name: str = Field(alias="Sample name")
    raw_reads: int = Field(alias="Raw reads")
    first_trim: int = Field(alias="First trim")
    second_trim: Optional[int] = Field(default=None, alias="Second trim")
    rRNA_mapped: Optional[int] = Field(default=None, alias="rRNA mapped")
    rRNA_free: Optional[int] = Field(default=None, alias="rRNA free")


class StarStats(BaseModel):
    """StarStats
    Pydantic model for STAR alignment statistics
    """

    input_reads: int = Field(default=0, alias="Input reads")
    # mapped reads
    mapped_total: int = Field(default=0, alias="Mapped: Total")
    mapped_unique_reads: Optional[int] = Field(
        default=None, alias="Mapped: Unique reads"
    )
    mapped_PCR_duplicate_reads: Optional[int] = Field(
        default=None, alias="Mapped: PCR duplicate reads"
    )
    mapped_uniquely_mapped_reads: int = Field(
        default=0, alias="Mapped: Uniquely mapped reads"
    )
    mapped_multimapped_reads: int = Field(default=0, alias="Mapped: Multimapped reads")
    # unmapped reads
    unmapped_total: Optional[int] = Field(
        default=None, alias="Unmapped: Total"
    )  # to accomodate deduplicated bams
    unmapped_mapped_to_too_many_loci: Optional[int] = Field(
        default=None, alias="Unmapped: mapped to too many loci"
    )
    unmapped_no_seed_windows: Optional[int] = Field(
        default=None, alias="Unmapped: no seed/windows"
    )
    unmapped_too_many_mismatches: Optional[int] = Field(
        default=None, alias="Unmapped: too many mismatches"
    )
    unmapped_too_short: Optional[int] = Field(default=None, alias="Unmapped: too short")
    unmapped_paired_end_mate: Optional[int] = Field(
        default=None, alias="Unmapped: paired-end mate"
    )
    model_config = ConfigDict(extra="allow")


class UmiStats(BaseModel):
    """UMIStats
    Pydantic model for  alignment statistics after UMI deduplication
    """

    umi_mapped_total: Optional[int] = Field(
        default=None, alias="Mapped: Total (UMI dedup)"
    )
    umi_mapped_uniquely_mapped_reads: Optional[int] = Field(
        default=None, alias="Mapped: Uniquely mapped reads (UMI dedup)"
    )
    umi_mapped_multimapped_reads: Optional[int] = Field(
        default=None, alias="Mapped: Multimapped reads (UMI dedup)"
    )
    model_config = ConfigDict(extra="allow")


class AllStats(UmiStats, StarStats, BaseStats):
    """AllStats
    Pydantic model combining BaseStats, Stats and UMI_Stats
    for comprehensive sample statistics from raw reads to alignment stats or deduplicated alignment stats.
    The inheritance order is important here output fields in proper order from raw to deduplicated stats
    """

    # @TODO: Add number validation
    model_config = ConfigDict(populate_by_name=True, extra="allow")
