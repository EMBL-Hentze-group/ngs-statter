import json
import logging
import re
from collections import defaultdict
from pathlib import Path
from random import sample

import numpy as np
from bokeh.layouts import column, row
from bokeh.models import (
    ColorPicker,
    ColumnDataSource,
    CustomJS,
    Legend,
    CheckboxGroup,
)
from bokeh.models.formatters import NumeralTickFormatter
from bokeh.palettes import Turbo256
from bokeh.plotting import figure, save

from jsonschema import validate

"""
Plot contamination vs unclassified read length distribuion
"""


class GanonReadLengthPlot:

    """
    Ganon contamination vs unclassified read length plotter
    """

    def __init__(
        self, json_folder: str, output_file: str, pattern: str = "json", nsamples=1
    ) -> None:
        self.json_folder = json_folder
        self.output_file = output_file
        self.pattern = pattern
        self.nsamples = nsamples
        self._read_lengths = set()
        self._read_length_dist = defaultdict(dict)
        self.samples = list()
        self._load_json()

    def _load_json(self) -> None:
        _json_dict = {}
        # schema to validate ganon read length json
        _json_schema = {
            "type": "object",
            "properties": {
                "contamination": {"type": "object"},
                "unclassified": {"type": "object"},
            },
            "required": ["contamination", "unclassified"],
        }
        for jf in sorted(Path(self.json_folder).glob(self.pattern)):
            sample = re.sub(r"\..*$", "", jf.stem)
            self.samples.append(sample)
            with open(jf, "r") as _jh:
                read_dat = json.load(_jh)
            validate(instance=read_dat, schema=_json_schema)
            # contamination read lengths
            self._read_lengths.update(
                list(map(lambda r: int(r), read_dat["contamination"]))
            )
            # unclassified read lengths
            self._read_lengths.update(
                list(map(lambda r: int(r), read_dat["unclassified"]))
            )
            _json_dict[sample] = read_dat
        for sample, read_len_dist in _json_dict.items():
            contamination = np.zeros((len(self._read_lengths)), dtype=np.uint32)
            unclassified = np.zeros_like(contamination)
            for i, rl in enumerate(sorted(self._read_lengths)):
                contamination[i] = read_len_dist["contamination"].get(str(rl), 0)
                unclassified[i] = read_len_dist["unclassified"].get(str(rl), 0)
            self._read_length_dist[sample] = (contamination, unclassified)

    def plot(self) -> None:
        rl_plot = figure(
            title=f"Read length distribution",
            x_axis_label="Read lengths",
            y_axis_label="Counts",
            width=1200,
            height=900,
        )
        legends, color_pickers, line_renders = list(), list(), list()
        colors = sample(Turbo256, len(self._read_length_dist))
        # line_counter =
        read_lengths = sorted(self._read_lengths)
        # 0: contamination, 1: unclassified
        checkbox_indices = {0: list(), 1: list()}
        line_counter = 0
        lines = list()
        sample_names = list()
        for i, sample_name in enumerate(sorted(self._read_length_dist)):
            i_rl = ColumnDataSource(
                data={
                    "rl": sorted(self._read_lengths),
                    "contamination": self._read_length_dist[sample_name][0],
                    "unclassified": self._read_length_dist[sample_name][1],
                }
            )
            # plot contamination read length distribution
            line_cont = rl_plot.line(
                x=read_lengths,
                y=self._read_length_dist[sample_name][0],
                line_width=2.5,
                line_alpha=0.8,
                line_dash="dashed",
                color=colors[i],
            )
            # plot unclassified read length distribution
            line_unc = rl_plot.line(
                x=read_lengths,
                y=self._read_length_dist[sample_name][1],
                line_width=2.5,
                line_alpha=0.8,
                color=colors[i],
            )
            sample_names.extend(
                [f"{sample_name} contamination", f"{sample_name} unclassified"]
            )
            line_renders.extend([line_cont, line_unc])
            # new line counter
            line_counter += 1
            # color picker
            picker = ColorPicker(title=sample_name, color=colors[i])
            picker.js_link("color", line_cont.glyph, "line_color")
            picker.js_link("color", line_unc.glyph, "line_color")
            color_pickers.append(picker)
            # legend
            legends.append((sample_name, [line_cont, line_unc]))
        # legend format
        # legend = Legend(items=legends, orientation="horizontal", ncols=self.nsamples)
        # legend.click_policy = "hide"
        # rl_plot.add_layout(legend, "below")
        # title
        rl_plot.title.align = "center"
        rl_plot.title.text_font_size = "16pt"
        rl_plot.title.text_font_style = "bold"
        # zoom limit
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
        print(sample_names)
        print(line_renders)
        # checkbox
        checkbox = CheckboxGroup(
            labels=sample_names, active=[i for i in range(len(sample_names))]
        )
        # custom js code
        js_code = CustomJS(
            args=dict(lines=line_renders, checkbox=checkbox),
            code="""
            for(var i=0; i<lines.length; i++){
                 lines[i].visible = checkbox.active.includes(i);
            }
            """,
        )
        checkbox.js_on_change("active", js_code)
        if Path(self.output_file).exists():
            logging.warning(f"Over writing {self.output_file}")
        # organize per sample color picker according to self.nsamples
        if self.nsamples > 1:
            col_picks = [
                row(*color_pickers[i : i + self.nsamples])
                for i in range(0, len(color_pickers), self.nsamples)
            ]
            save(
                row(rl_plot, column(column(*col_picks), checkbox)),
                filename=self.output_file,
                title="Gene type read length distribution",
            )
        else:
            save(
                row(rl_plot, column(column(*color_pickers), checkbox)),
                filename=self.output_file,
                title="Gene type read length distribution",
            )


grlp = GanonReadLengthPlot(
    "/workspaces/mad_statter/test_data/",
    "/workspaces/mad_statter/test_data/test_output.html",
    "*read_length.json",
)

grlp.plot()
