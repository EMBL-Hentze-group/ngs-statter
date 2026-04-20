import logging
import sys
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)


class YamlCheck:
    """
    Check sanity of the supplied yaml file
    """

    def __init__(self):
        # each file must have the following keys
        # name: name to call the sample
        # group: group to which the sample belongs
        self.groups = set(["name", "group"])

    def check_sanity(self, fin):
        sanity_map = {}
        with open(fin, "r") as yh:
            ymap = yaml.safe_load(yh)
        for k, v in ymap.items():
            pk = Path(k)
            if (not pk.exists()) or (not pk.is_file()):
                logger.warning(
                    "{}: either this path is not a file or it does not exists! skipping...".format(
                        k
                    )
                )
                continue
            diffk = self.groups - set(v.keys())
            if len(diffk) > 0:
                logger.warning(
                    "Cannot find attribute(s) {} for {}! skipping....".format(
                        ", ".join(diffk), k
                    )
                )
            logger.info(
                "File {},  sample name: {}, group: {} ".format(k, v["name"], v["group"])
            )
            sanity_map[k] = v
        if len(sanity_map) <= 1:
            error = """
            Only 1 or no file(s) has been successfully parsed from {}!
            Make sure that the input yaml file follows the format:

            /path/to/input_1_file.ext
                name: IP1 # a suitable name for the sample
                group: IP # group to which sample belongs
            /path/to/input_2_file.ext
                name: SMI1 # a suitable name for the sample
                group: SMI
            ....

            """.format(
                fin
            )
            raise RuntimeError(error)
        return sanity_map

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
