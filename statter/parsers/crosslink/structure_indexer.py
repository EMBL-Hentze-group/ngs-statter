import logging
from venv import logger

import numpy as np
from rtree import index

from statter.parsers.crosslink.sites_to_bed import SitesToBed

logger = logging.getLogger(__name__)


class StructureIndexer:
    """
    Parse given secondary structure/primary motif file
    extend 5' and 3' regions, return as rtree index
    """

    def __init__(self, bed, l=100, r=100, unstranded=False) -> None:
        self.bed = bed
        self.l = l
        self.r = r
        self.unstranded = unstranded
        self.index = {}

    def _get_slop(self):
        """
        return appropriate function to slop input region
        """

        def do_nothing(begin, end, strand):
            # if self.l and self.r are 0
            return (begin, end)

        def unstranded(begin, end, strand):
            # if unstranded = True
            return (begin - self.l, end + self.r)

        def stranded(begin, end, strand):
            # unstranded = False and either
            if strand == "-":
                newBegin = max(0, begin - self.r)
                return (newBegin, end + self.l)
            else:
                return (begin - self.l, end + self.r)

        if (self.l == 0) and (self.r == 0):
            logger.debug("No extension l = 0 and r = 0")
            return do_nothing
        elif self.unstranded:
            logger.debug("extend unstranded, l = {} and r = {}".format(self.l, self.r))
            return unstranded
        else:
            logger.debug("extend stranded, l = {} and r = {}".format(self.l, self.r))
            return stranded

    def extend_index_regions(self):
        freader, mode = SitesToBed._get_parser(self.bed)
        distlist, extended_distlist = list(), list()
        slopper = self._get_slop()
        #  slopped_out = str(self.tmpdir/"struct_l{}_r{}.bed".format(self.l, self.r))
        iids = 0
        with freader(self.bed, mode) as fh:  # , open(slopped_out,'w') as oh :
            for ix, f in enumerate(fh):  # in case name is not unique per entry
                fdat = f.strip().split("\t")
                try:
                    begin = int(fdat[1])
                except ValueError as v:
                    logger.warning(str(v) + " skipping....")
                    continue
                end = int(fdat[2])
                distlist.append(end - begin)
                new_begin, new_end = slopper(begin, end, fdat[-1])
                extended_distlist.append(new_end - new_begin)
                if fdat[5] == "+":
                    strand = (0, 0)
                elif fdat[5] == "-":
                    strand = (1, 1)
                else:
                    strand = (0, 1)
                dat_map = {
                    "name": "{}_{}".format(fdat[3], ix),
                    "strand": fdat[5],
                    "begin": begin,
                    "end": end,
                }
                logger.debug(dat_map)
                try:
                    self.index[fdat[0]].insert(
                        iids, (new_begin, strand[0], new_end, strand[1]), dat_map
                    )
                except KeyError:
                    self.index[fdat[0]] = index.Index()
                    self.index[fdat[0]].insert(
                        iids, (new_begin, strand[0], new_end, strand[1]), dat_map
                    )
                iids += 1
                # oh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(fdat[0], new_begin, new_end, fdat[3], fdat[4],+ fdat[5]))
        logger.info("{} parsed {} lines, indexed {} regions".format(self.bed, ix, iids))
        return (
            self.index,
            np.array(distlist, dtype=np.float32),
            np.array(extended_distlist, dtype=np.float32),
        )
