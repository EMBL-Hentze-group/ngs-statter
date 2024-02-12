import json
import logging
import re
from pathlib import Path

import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColorPicker, Legend
from bokeh.models.formatters import NumeralTickFormatter
from bokeh.plotting import figure, save
from bokeh.palettes import Turbo256
import random

"""
plot read length distribution after adapter trimming
@TODO: add vertical line for min. length
"""


class ReadLengthPlot:
    def __init__(
        self,
        json_folder: str,
        output_file: str,
        min_read_length: int = 15,
        pattern: str = "*.json",
        title: str = "Read length distribution after adapter trimming",
        nsamples: int = 1,
    ) -> None:
        self.json_folder = json_folder
        self.pattern = pattern
        self.min_read_length = min_read_length
        self.output_file = output_file
        self.read_lens = {}
        self._get_read_lengths()
        self.title = title
        self.nsamples = nsamples

    def _get_read_lengths(self) -> None:
        """
        Helper function
        Read json formatted read lengths
        """
        for jp in sorted(Path(self.json_folder).glob(self.pattern)):
            with open(jp, "r") as jh:
                self.read_lens[re.sub("\..*$", "", jp.stem)] = json.load(jh)
        if len(self.read_lens) == 0:
            raise ValueError(
                f"Cannot find json formatted read length files in folder {self.json_folder}!"
            )

    def plot(self) -> None:
        read_lengths = set()
        for s, dat in self.read_lens.items():
            read_lengths.update(dat.keys())
        read_lengths = sorted([int(r) for r in read_lengths])
        rl_plot = figure(
            title=self.title,
            x_axis_label="Read lengths",
            y_axis_label="Counts",
            width=1200,
            height=900,
        )
        if self.min_read_length > 0:
            rl_plot.vspan(
                x=self.min_read_length,
                line_width=1.5,
                line_dash="dashed",
                line_color="#555753",
            )
        legends, col_pickers = list(), list()
        colors = random.sample(Turbo256, len(self.read_lens))
        for c, sample in enumerate(sorted(self.read_lens.keys())):
            sample_name = re.sub(r"(\_{1,}|\-{1,}|\.{1,})", " ", sample)
            counts = np.zeros(len(read_lengths), dtype=np.int32)
            for i, rl in enumerate(read_lengths):
                try:
                    counts[i] = self.read_lens[sample][str(rl)]
                except KeyError:
                    pass
            line = rl_plot.line(
                x=read_lengths,
                y=counts,
                line_width=2.5,
                line_alpha=0.8,
                color=colors[c],
            )
            legends.append((sample_name, [line]))
            picker = ColorPicker(title=sample_name, color=colors[c])
            picker.js_link("color", line.glyph, "line_color")
            col_pickers.append(picker)
        legend = Legend(items=legends, orientation="horizontal", ncols=self.nsamples)
        legend.click_policy = "hide"
        rl_plot.add_layout(legend, "below")
        # title
        rl_plot.title.align = "center"
        rl_plot.title.text_font_size = "16pt"
        rl_plot.title.text_font_style = "bold"
        # limit zoom limit
        rl_plot.x_range.min_interval = 1
        # axis labels
        # X
        rl_plot.xaxis.major_label_text_font_size = "12pt"
        rl_plot.xaxis.axis_label_text_font_size = "14pt"
        rl_plot.xaxis.axis_label_text_font_style = "bold"
        # Y
        rl_plot.yaxis.major_label_text_font_size = "12pt"
        rl_plot.yaxis.axis_label_text_font_size = "14pt"
        rl_plot.yaxis.axis_label_text_font_style = "bold"
        # add comma to counts on Y axis
        rl_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
        # organize color picker for samples
        if Path(self.output_file).exists():
            logging.warning(f"Over writing {self.output_file}")
        if self.nsamples > 1:
            col_picks = [
                row(*col_pickers[i : i + self.nsamples])
                for i in range(0, len(col_pickers), self.nsamples)
            ]
            save(
                row(rl_plot, column(*col_picks)),
                filename=self.output_file,
                title="Read length distribution",
            )
        else:
            save(
                row(rl_plot, column(*col_pickers)),
                filename=self.output_file,
                title="Read length distribution",
            )
