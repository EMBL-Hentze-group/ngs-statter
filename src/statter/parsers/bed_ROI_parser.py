import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


class RegionReader:
    """
    Parse given secondary structure/primary motif file
    extend 5' and 3' regions, return as rtree index
    """

    def __init__(
        self,
        bed: str | Path,
        l: int = 100,
        r: int = 100,
        unstranded: bool = False,
        most_5prime: bool = False,
    ) -> None:
        self._bed = self._read_bed(bed, l, r, unstranded)
        self.index = {}
        if most_5prime:
            self._bed = self._5prime_most()

    def _read_bed(
        self, bedf: str | Path, l: int, r: int, unstranded: bool
    ) -> pd.DataFrame:
        """_read_bed read BED file and extend regions according to l and r parameters

        Args:
            bedf: (str) path to BED file
            l: (int) number of bases to extend at the 5' end
            r: (int) number of bases to extend at the 3' end
            unstranded: (bool) whether to ignore strand information

        Raises:
            ValueError: if the BED file has fewer than 6 columns

        Returns:
            pd.DataFrame: DataFrame containing the extended BED regions
        """
        names = ["chrom", "start", "end", "name", "score", "strand"]
        bed = pd.read_csv(bedf, sep="\t", header=None, names=names)
        if bed.shape[1] < 6:
            raise ValueError(f"BED file {bedf} must have at least 6 columns!")
        if unstranded:
            bed["start_extend"] = (bed["start"] - l).clip(lower=0)
            bed["end_extend"] = bed["end"] + r
        else:
            bed["start_extend"] = bed.apply(
                lambda row: (
                    max(0, row["start"] - r)
                    if row["strand"] == "-"
                    else max(0, row["start"] - l)
                ),
                axis=1,
            )
            bed["end_extend"] = bed.apply(
                lambda row: (
                    row["end"] + l if row["strand"] == "-" else row["end"] + r
                ),
                axis=1,
            )
        return bed

    def _5prime_most(self) -> pd.DataFrame:
        """_5prime_most Helper function
        If the regions after extension overlaps, pick the most 5' region relative to strand out of the overlapping regions
        Returns:
            pd.DataFrame: DataFrame containing the selected 5' most regions
        """
        _5pmost: list[pd.DataFrame] = list()
        chrom_groups = self._bed.groupby(["chrom", "strand"])
        for (chrom, strand), gdf in chrom_groups:
            sort_col: str = "start_extend"
            if strand == "-":
                sort_col = "end_extend"
            _5pdf = _pick_5prime_most(
                gdf,
                strand=strand,
                sort_col=sort_col,
                extended_start_col="start_extend",
                extended_end_col="end_extend",
            )
            _5pmost.append(_5pdf)
        return pd.concat(_5pmost)

    @property
    def regions(self) -> pd.DataFrame:
        """
        return the extended regions as a DataFrame
        """
        return self._bed

    @property
    def max_length(self):
        """
        return maximum length of the regions before extension
        """
        return (self._bed["end"] - self._bed["start"]).max()

    # def extend_index_regions(self):
    #     freader, mode = SitesToBed._get_parser(self.bed)
    #     distlist, extended_distlist = list(), list()
    #     slopper = self._get_slop()
    #     #  slopped_out = str(self.tmpdir/"struct_l{}_r{}.bed".format(self.l, self.r))
    #     iids = 0
    #     with freader(self.bed, mode) as fh:  # , open(slopped_out,'w') as oh :
    #         for ix, f in enumerate(fh):  # in case name is not unique per entry
    #             fdat = f.strip().split("\t")
    #             try:
    #                 begin = int(fdat[1])
    #             except ValueError as v:
    #                 logger.warning(str(v) + " skipping....")
    #                 continue
    #             end = int(fdat[2])
    #             distlist.append(end - begin)
    #             new_begin, new_end = slopper(begin, end, fdat[-1])
    #             extended_distlist.append(new_end - new_begin)
    #             if fdat[5] == "+":
    #                 strand = (0, 0)
    #             elif fdat[5] == "-":
    #                 strand = (1, 1)
    #             else:
    #                 strand = (0, 1)
    #             dat_map = {
    #                 "name": "{}_{}".format(fdat[3], ix),
    #                 "strand": fdat[5],
    #                 "begin": begin,
    #                 "end": end,
    #             }
    #             logger.debug(dat_map)
    #             try:
    #                 self.index[fdat[0]].insert(
    #                     iids, (new_begin, strand[0], new_end, strand[1]), dat_map
    #                 )
    #             except KeyError:
    #                 self.index[fdat[0]] = index.Index()
    #                 self.index[fdat[0]].insert(
    #                     iids, (new_begin, strand[0], new_end, strand[1]), dat_map
    #                 )
    #             iids += 1
    #             # oh.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(fdat[0], new_begin, new_end, fdat[3], fdat[4],+ fdat[5]))
    #     logger.info("{} parsed {} lines, indexed {} regions".format(self.bed, ix, iids))
    #     return (
    #         self.index,
    #         np.array(distlist, dtype=np.float32),
    #         np.array(extended_distlist, dtype=np.float32),
    #     )


def _pick_5prime_most(
    gdf: pd.DataFrame,
    strand: str = "+",
    sort_col: str = "start_extend",
    extended_start_col: str = "start_extend",
    extended_end_col: str = "end_extend",
) -> pd.DataFrame:
    """_pick_5prime_most If there are overlapping regions, pick the most 5' region relative to strand
    Args:
        gdf: (pd.DataFrame) DataFrame containing the regions to be filtered
        strand: (str) strand information, either "+" or "-". Defaults to "+".
        sort_col: (str) column name to sort by to determine 5' most region.
        extended_start_col: (str) column name for extended start positions. Defaults to "start_extend".
        extended_end_col: (str) column name for extended end positions. Defaults to "end_extend".

    Returns:
        pd.DataFrame: DataFrame containing the selected 5' most regions.
    """

    def plus_strand_check(current, rnext, extended_start_col, extended_end_col):
        return current[extended_end_col] <= rnext[extended_start_col]

    def minus_strand_check(current, rnext, extended_start_col, extended_end_col):
        return not (current[extended_start_col] <= rnext[extended_end_col])

    if strand == "-":
        gdf = gdf.sort_values(sort_col, ascending=False)
        checkfn = minus_strand_check
    else:
        gdf = gdf.sort_values(sort_col, ascending=True)
        checkfn = plus_strand_check
    if gdf.shape[0] == 1:
        return gdf
    else:
        selected = list()
        idf = gdf.iterrows()
        try:
            _, current = next(idf)
            for _, rnext in idf:
                if checkfn(current, rnext, extended_start_col, extended_end_col):
                    selected.append(current)
                    current = rnext
            selected.append(current)
        except StopIteration:
            pass
        _df = pd.DataFrame(selected)
        return _df.sort_values(extended_start_col, ascending=True)
