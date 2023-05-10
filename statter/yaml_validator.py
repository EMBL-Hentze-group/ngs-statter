import yaml
from typing import Set

"""
Module(s) to validate input yaml file(s)
"""


class BaseYaml:
    """
    Base yaml file follows the following format:
    sample name:
        path: /path/to/file
        group: sample group name
    """

    def __init__(self, yaml_file: str) -> None:
        self._required_keys: Set[str] = set(["path", "group"])
        with open(yaml_file, "r") as _yh:
            self._yaml_data = yaml.load(_yh, yaml.FullLoader)

    def _base_validator(self) -> None:
        """
        Helper function
        Validate base yaml data
        """
        for label, dat in self._yaml_data.items():
            missing_keys = self._required_keys - set(dat.keys())
            if len(missing_keys) > 0:
                msg = f"Cannot find all required keys for sample {label}. Required keys: {', '.join(self._required_keys)}. Missing key(s): {', '.join(missing_keys)}"
                raise RuntimeError(msg)
            if not isinstance(dat["path"], str):
                msg = f"Sample {label}, 'path' value MUST BE a string '/path/to/file', but found {type(dat['path'].__name__)}"
                raise RuntimeError(msg)
            if not isinstance(dat["group"], str):
                msg = f"Sample {label}, 'group' value MUST BE a string name of the group that this sample belongs to, but found {type(dat['group'].__name__)}"
                raise RuntimeError(msg)
