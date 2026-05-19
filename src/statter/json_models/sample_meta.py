import logging
from pathlib import Path
from typing import Optional

from matplotlib.colors import is_color_like
from pydantic import BaseModel, ConfigDict, Field, field_validator

logger = logging.getLogger(__name__)


class SampleMeta(BaseModel):
    """SampleMeta
    Pydantic model for sample metadata
    """

    file: str = Field(alias="File")
    sample: str = Field(alias="Sample")
    group: str = Field(alias="Group")
    color: Optional[str] = Field(default=None, alias="Color")

    model_config = ConfigDict(
        extra="allow",
        populate_by_name=True,
    )

    @field_validator("file")
    def validate_file(cls, file: str) -> str:
        fp = Path(file)
        if not fp.is_file():
            raise ValueError(f"File '{file}' must be a valid file path.")
        elif not fp.exists():
            raise ValueError(f"File '{file}' does not exist.")
        return str(fp.resolve())

    @field_validator("color")
    def validate_color(cls, color: Optional[str]) -> Optional[str]:
        if (color is not None) and (not is_color_like(color)):
            logger.warning(
                f"Color '{color}' is not a valid color. Setting color to None, all colors will be generated randomly."
            )
            return None
        return color

    @classmethod
    def required_fields(cls) -> set[str]:
        cols: set[str] = set()
        for k, field in cls.model_fields.items():
            if field.is_required():
                cols.add(k)
        return cols
