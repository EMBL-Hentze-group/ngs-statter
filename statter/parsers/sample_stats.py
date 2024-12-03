import logging
import json
from typing import Dict, List
from collections import OrderedDict
import re
from pathlib import Path
import pandas as pd
from pysam import Union

logger = logging.getLogger(__name__)
class SampleStats:
    def __init__(self, first_trim:str, second_trim:str, rRNA_free:str, rRNA_mapped:str, align_stats:str) -> None:
        self.first_trim = first_trim
        self.second_trim = second_trim
        self.rRNA_free = rRNA_free
        self.rRNA_mapped = rRNA_mapped
        self.align_stats = align_stats
        
        self._sample_stats:Dict = OrderedDict()
        # collect stats
        # TODO: json schema
    
    @staticmethod
    def _load_json(fname: str) -> Dict:
        """load json

        Helper function, load json file.

        Args:
            fname (str): Json file name

        Returns:
            Dict: contents of the json file
        """
        with open(fname, "r") as jh:
            dat = json.load(jh)
        return dat
    
    def _add_trimming_stats(self) -> None:
        """_add_trimming_stats Helper function
        add trimming stats to the sample stats dictionary
        """
        first_trim = self._load_json(self.first_trim)
        second_trim = self._load_json(self.second_trim)
        omit_keys = set(["filtering_result"])
        self._sample_stats["Raw reads"] = first_trim["summary"]["before_filtering"]["total_reads"]
        for comb in (first_trim["filtering_result"],second_trim["filtering_result"]):
            for k,v in comb.items():
                if k in omit_keys:
                    continue
                stat_name = re.sub(r"\_{1,}"," ",k)
                try:
                    self._sample_stats[stat_name]+=v
                except KeyError:
                    self._sample_stats[stat_name]=v
        self._sample_stats["Trimmed reads"] = second_trim["summary"]["after_filtering"]["total_reads"]
    
    def _add_rRNA_stats(self) -> None:
        """ _add_rRNA_stats Helper function
        add rRNA stats to the sample stats dictionary.
        """
        rRNA_mapped = sum(self._load_json(self.rRNA_mapped).values())
        rRNA_free = sum(self._load_json(self.rRNA_free).values())
        if self._sample_stats["Trimmed reads"]!= rRNA_mapped + rRNA_free:
            raise RuntimeError("The number of trimmed reads does not match the number of mapped rRNA + unmapped rRNA reads!")
        self._sample_stats["rRNA mapped"] = rRNA_mapped
        self._sample_stats["rRNA free"] = rRNA_free
        
    def  _add_alignment_stats(self) -> None:
        """_add_alignment_stats Helper function
        add alignment stats to the sample stats dictionary.
        """
        align_stats = self._load_json(self.align_stats)
        if align_stats["Input reads"] != self._sample_stats["rRNA free"]:
            raise RuntimeError("The number of rRNA reads does not match the number of genome mapped reads!")
        for k,v in align_stats.items():
            if re.match(r"^.*\%$",k):
                continue
            self._sample_stats[k]=v
        # self._sample_stats[""] = self._load_json(self.rRNA_free)["summary"][0]["total_reads"]
        # self._sample_stats["Mapped reads"] = self._load_json(self.rRNA_mapped)["summary"][0]["total_reads"]
    
    def collect_stats(self,out:str) -> None:
        """collect_stats 
        Collect all stats and write to a json file
        Args:
            out: output file name
        """
        if Path(out).exists():
            logger.warning(f"Re-writing file {out}")
        self._add_alignment_stats()
        self._add_rRNA_stats()
        self._add_alignment_stats()
        with open("out","w") as _jh:
            json.dump(self._sample_stats, _jh)

def gather_sample_stats(stats_dir:Union[str,Path],out:str,suffix="sample_stats.json",) -> None:
    """gather_sample_stats 
    Gather all sample stats with given prefix and write them into a single csv file

    Args:
        stats_dir: directory with sample stats json files
        out: output csv file name
        suffix: sample stats json suffix. Defaults to "sample_stats.json".
    """
    samples:List[Path] = sorted(Path(stats_dir).glob(f"{suffix}"))
    if len(samples) == 0:
        raise RuntimeError(f"Cannot find sample stats file with suffix {suffix} in directory {str(stats_dir)}")
    stats,names = [],[]
    for s in samples:
        with open(s,"r") as _jh:
            stats.append(json.load(_jh))
        names.append(re.sub(r"\..*$","",s.stem))
    if Path(out).exists():
        logger.warning(f"Re-writing file {out}")
    all_stats = pd.DataFrame(stats,index=names)
    all_stats.to_csv(sep="\t",index=True)