import logging
import tempfile
from pathlib import Path
from shutil import rmtree

from statter.parsers.crosslink.overlap_finder import FindOverlaps
from statter.parsers.crosslink.sites_to_bed import SitesToBed
from statter.parsers.crosslink.structure_indexer import StructureIndexer
from statter.parsers.crosslink.yaml_parser import YamlCheck

logger = logging.getLogger(__name__)


class Plotter:
    """
    Class to wrap everything together
    """

    def __init__(self, yaml, region_bed, tmpdir=None) -> None:
        self.yaml = yaml
        self.region_bed = region_bed
        if tmpdir is None:
            self.tmpdir = Path(tempfile.gettempdir())
            logger.debug("Using system tmp folder: {}".format(self.tmpdir))
        else:
            tempfile.tempdir = tmpdir
            self.tmpdir = Path(tempfile.tempdir)
            logger.debug("Using user provided tmp folder: {}".format(str(self.tmpdir)))
        self._file_map = {}  # formatted and aggregated file map
        self._index = {}  # indexed regions
        self.l = 0  # 5' extension
        self._orig_leng = None  # np.array of original lengths
        self._extended_leng = None  # np.array of extended lengths

    def __enter__(self):
        yc = YamlCheck()
        self._in_yml = yc.check_sanity(self.yaml)
        if not self.tmpdir.exists():
            self.tmpdir.mkdir()
        self._tmpdir = self.tmpdir / next(tempfile._get_candidate_names())  # type: ignore
        logger.info("Using {} for temp files".format(str(self._tmpdir)))
        return self

    def __exit__(self, except_type, except_val, except_traceback):
        rmtree(self._tmpdir)
        if except_type:
            logger.exception(except_val)

    def aggregate_sites(self, norm=None):
        stb = SitesToBed(self._tmpdir, self._in_yml, norm)
        self._file_map = stb.aggregate()

    def index_region(self, l=100, r=100):
        self.l = l
        si = StructureIndexer(self.region_bed, l=l, r=r)
        self._index, self._orig_leng, self._extended_leng = si.extend_index_regions()

    def plot(
        self,
        out_file,
        smoothing_window=5,
        xlabel="Relative crosslink positions",
        ylabel="crosslink count",
        width=30,
        height=27,
        errorbar="sd",
    ):
        if errorbar == "None":
            errorbar = None
        fo = FindOverlaps(
            self._index, self._orig_leng, self._extended_leng, self._file_map, self.l
        )
        fo.plot(out_file, smoothing_window, xlabel, ylabel, width, height, errorbar)
