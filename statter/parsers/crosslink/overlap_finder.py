import logging
from operator import itemgetter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)


class FindOverlaps:
    """
    Find crosslink site counts that overlap given region of interest
    count and aggregate crosslink site count data and finally plot
    """

    def __init__(
        self, index, original_lengths, region_lengths, file_map, l=100
    ) -> None:
        self.index = index
        self.files = file_map
        self.orig_len = np.array(original_lengths)
        self.reg_len = np.array(region_lengths)
        self.max_array = int(np.max(self.reg_len))
        self.l = l
        logger.info("Maximum number of data points is {}".format(self.max_array))
        self.dat = list()

    def _get_smoother(self, window_size=1):

        def do_nothing(counts):
            return counts

        def smoother(counts):
            """
            smooth using a moving average
            """
            if window_size > len(counts):
                logger.warning(
                    "Normalization window ({}) is larger than the number of counts (size: {}) !. Using non normalized counts ".format(
                        window_size, len(counts)
                    )
                )
                return counts
            elif window_size == 2:
                kernel = np.array([0.5, 0.5])
            else:
                # see https://stackoverflow.com/questions/17190649/how-to-obtain-a-gaussian-filter-in-python
                k = np.arange(-window_size // 2 + 1, window_size // 2 + 1, 1)
                knorm = np.exp(-(k**2) / (2.0**2))
                kernel = knorm / knorm.sum()
            return np.convolve(counts, kernel, mode="same")

        if window_size <= 1:
            return do_nothing
        else:
            logger.info(
                "Using gaussian smoother, smoothing using {} adjacent positions".format(
                    window_size
                )
            )
            return smoother

    def get_counts(self, smoothing_window=1):
        """
        Compute counts per individual nucleotide position
        """
        smoother = self._get_smoother(window_size=smoothing_window)
        for f, d in self.files.items():
            pos_dict = {}
            skiplist = set()
            with open(f, "r") as fh:
                logger.info(
                    "Working on sample: {} group: {}".format(d["name"], d["group"])
                )
                logger.debug("Temp file {}".format(f))
                for f in fh:
                    fdat = f.strip().split("\t")
                    if fdat[0] not in self.index:
                        if fdat[0] not in skiplist:
                            logger.info(
                                "Given bed has no structure in {}, skipping....".format(
                                    fdat[0]
                                )
                            )
                            skiplist.add(fdat[0])
                        continue
                    try:
                        begin = int(fdat[1])
                    except ValueError as v:
                        logger.warning(str(v) + " skipping...")
                        continue
                    end = int(fdat[2])
                    score = float(fdat[4])
                    if fdat[5] == "+":
                        strand = (0, 0)
                    elif fdat[5] == "-":
                        strand = (1, 1)
                    else:
                        strand = (0, 1)
                    for ix in self.index[fdat[0]].intersection(
                        (begin, strand[0], end, strand[1]), objects=True
                    ):
                        if ix.bbox[2] == begin:
                            # 0th position base of the crosslink and last base of region overlap,
                            # skip
                            continue
                        if fdat[5] == "-":
                            pos = int(ix.bbox[2] - end)
                        else:
                            pos = int(end - ix.bbox[0])
                        if pos >= self.max_array:
                            logger.info(
                                "Skipping position {} in {} {}".format(
                                    end, d["name"], d["group"]
                                )
                            )
                            continue
                        logger.debug("{} {} {}".format(ix.object["name"], pos, score))
                        try:
                            pos_dict[ix.object["name"]][pos] = score
                        except KeyError:
                            pos_dict[ix.object["name"]] = [0] * self.max_array
                            pos_dict[ix.object["name"]][pos] = score
            regions = sorted(pos_dict.keys())
            for r in regions:
                for i, score in enumerate(smoother(pos_dict[r])):
                    self.dat.append([i, score, d["name"], d["group"]])

    def _plot(
        self,
        file,
        xlabel="Relative position",
        ylabel="crosslink count",
        width=30,
        height=27,
        errorbar="sd",
    ):
        """
        plot data using seaborn
        """
        sns.reset_orig()
        # sns.set_style("white")
        df = pd.DataFrame(
            self.dat, columns=["relative_position", "counts", "sample", "group"]
        )
        # df.to_pickle(file+'.pkl')
        logger.debug("Default figure size: {}".format(plt.rcParams["figure.figsize"]))
        cm = 1 / 2.54  # centimeters in inches
        plt.rcParams["figure.figsize"] = (width * cm, height * cm)
        sns.set(rc={"figure.figsize": (width * cm, height * cm)})
        logger.debug("New figure size: {}".format(plt.rcParams["figure.figsize"]))
        sns.set_style("ticks")
        fig, ax = plt.subplots()
        sg = sns.lineplot(
            data=df,
            x="relative_position",
            y="counts",
            hue="group",
            style="sample",
            dashes=False,
            errorbar=errorbar,
        )
        sg.axvline(self.l, 0.02, 0.98, color="gray", ls="--")
        struct_90 = self.l + np.ceil(np.quantile(self.orig_len, 0.90))
        sg.axvline(struct_90, 0.02, 0.98, color="gray", ls="--")
        # format x axis position labels and add major ticks at every 10th position within the structure
        tick_pos = set(sg.get_xticks()[1:-1])
        tick_pos.update(np.arange(self.l, struct_90, 10))
        tick_pos = sorted(tick_pos)
        # remove trailing zeors in labels
        tick_labels = list(map(lambda tp: int(tp - self.l), tick_pos))
        sg.set_xticks(tick_pos)
        sg.set_xticklabels(tick_labels)
        sg.set_xlabel(xlabel)
        sg.set_ylabel(ylabel + " ({} sequences)".format(self.orig_len.shape[0]))
        # add minor ticks at every 5th position within the structure
        ax.set_xticks(np.arange(self.l + 5, struct_90 - 5, 5), minor=True)
        # remove "sample" from legend
        h, l = sg.get_legend_handles_labels()
        uniq_grps = len(df["group"].unique()) + 1
        # remove "group" heading from legend and use only the "group" labelling
        sg.legend(h[1:uniq_grps], l[1:uniq_grps])
        sns.despine()
        sfig = sg.get_figure()
        # if file.endswith('pdf') or file.endswith('svg'):
        #     sfig.savefig(file,bbox_inches='tight')
        # else:
        sfig.savefig(file, bbox_inches="tight", dpi=300)
        out_csv = str(Path(file).with_suffix(".csv"))
        self._save_csv(out_csv)

    def plot(
        self,
        file,
        smoothing_window=5,
        xlabel="Relative position",
        ylabel="crosslink count",
        width=30,
        height=27,
        errorbar="sd",
    ):
        """
        wrapper to call get_counts and plot results
        """
        self.get_counts(smoothing_window)
        self._plot(
            file=file,
            xlabel=xlabel,
            ylabel=ylabel,
            width=width,
            height=height,
            errorbar=errorbar,
        )

    def _save_csv(self, csv_out):
        """
        Helper function
        Write the corresponding  data to csv file
        """
        dat_map = {}
        for d in self.dat:
            try:
                dat_map[(d[0], d[2], d[3])].append(d[1])
            except KeyError:
                dat_map[(d[0], d[2], d[3])] = [d[1]]
        if Path(csv_out).exists():
            logger.warning(f"Re writing {csv_out}")
        with open(csv_out, "w") as fh:
            fh.write("Rel.position\tSample_name\tSample_type\tmean_count\n")
            for k in sorted(dat_map.keys(), key=itemgetter(0, 1)):
                dat = dat_map[k]
                mdat = np.mean(dat)
                # if np.sum(dat) == 0:
                #     ci = (0,0)
                # else:
                #     ci = st.t.interval(confidence=0.95, df=len(dat)-1, loc=mdat, scale=st.sem(dat))
                fh.write(f"{k[0]-self.l}\t{k[1]}\t{k[2]}\t{mdat}\n")
