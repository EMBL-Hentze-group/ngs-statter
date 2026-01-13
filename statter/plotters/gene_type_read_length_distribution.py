import json
import logging
import re
from collections import defaultdict
from itertools import chain
from pathlib import Path
from random import sample
from typing import Dict, List, Optional

from bokeh.layouts import column, row
from bokeh.models import ColorPicker, ColumnDataSource, CustomJS, Legend, Select
from bokeh.models.formatters import NumeralTickFormatter
from bokeh.palettes import Turbo256
from bokeh.plotting import figure, save

"""
Plot aligned read length distribution per gene type
"""


class GeneTypePlot:
    def __init__(
        self,
        output_file: str,
        json_files: Optional[List[str]] = None,
        json_folder: Optional[str] = None,
        pattern: str = "*.json",
        nsamples=1,
    ) -> None:
        self.json_folder = json_folder
        self.pattern = pattern
        self.output_file = output_file
        self.nsamples = nsamples
        self._gene_type_read_lengths = {}  # read length per gene type
        self._libsizes = {}  # library size for each sample
        self.source_data = {}
        self.samples = list()
        self._jsons: list[Path] = []
        if len(json_files) > 0:  # type: ignore
            self._jsons = sorted([Path(jf) for jf in json_files])  # type: ignore
        elif json_folder is not None:
            self._jsons = sorted(Path(json_folder).glob(pattern))
        else:
            raise RuntimeError("Either json_folder or json_files must be provided")
        self._get_data()

    def _get_data(self) -> None:
        """
        Helper function
        parse and format data
        """
        _tmp_data_dict = {}
        for jf in sorted(self._jsons):
            jf_name = re.sub(r"\..*$", "", jf.stem)
            self.samples.append(jf_name)
            _tmp_data_dict[jf_name] = self._load_format_convert(jf)

        for gene_type in sorted(self._gene_type_read_lengths.keys()):
            read_length_cds = list()
            read_lengths = sorted(self._gene_type_read_lengths[gene_type])
            nskipped = 0
            sample_lib_sizes = [0] * len(self.samples)
            for i, sample in enumerate(sorted(self.samples)):
                counts = [0] * len(read_lengths)
                if not gene_type in _tmp_data_dict[sample]:
                    logging.warning(
                        f"Cannot find gene type {gene_type} in sample {sample}, filling with 0s"
                    )
                    read_length_cds.append((read_lengths, counts))
                    nskipped += 1
                    continue
                for k, rl in enumerate(read_lengths):
                    if not rl in _tmp_data_dict[sample][gene_type]:
                        continue
                    counts[k] = _tmp_data_dict[sample][gene_type][rl]

                read_length_cds.append((read_lengths, counts))
                sample_lib_sizes[i] = sum(counts)
            if nskipped == len(self.samples):
                # do not add gene types if its fully populates with 0s
                continue
            self._libsizes[gene_type] = sample_lib_sizes
            self.source_data[gene_type] = read_length_cds

    def _load_format_convert(self, json_file: Path) -> Dict[str, Dict[int, int]]:
        """
        Helper function
        add "all" gene type
        convert read lengths from string to int
        """
        with open(json_file, "r") as jh:
            json_dat = json.load(jh)
        json_dict = defaultdict(dict)
        for gene_type, dat in json_dat.items():
            read_lengths = set()
            for rls, counts in dat.items():
                json_dict[gene_type][int(rls)] = counts
                try:
                    json_dict["all"][int(rls)] += counts
                except KeyError:
                    json_dict["all"][int(rls)] = counts
                read_lengths.add(int(rls))
            try:
                self._gene_type_read_lengths["all"].update(read_lengths)
            except KeyError:
                self._gene_type_read_lengths["all"] = read_lengths
            try:
                self._gene_type_read_lengths[gene_type].update(read_lengths)
            except KeyError:
                self._gene_type_read_lengths[gene_type] = read_lengths
        return json_dict

    def _group_n_sort_gene_types(self) -> List[str]:
        """
        Helper function
        group self._gene_type_read_lengths.keys() according to the endings ("_tRNA" or "_gene" or "_pseudogene" or None)
        """
        tRNA_regex = re.compile(r"^.*(\s{1,}|\_{1,})tRNA$", re.IGNORECASE)
        gene_regex = re.compile(r"^.*(\s{1,}|\_{1,})gene$", re.IGNORECASE)
        pseudogene_regex = re.compile(r"^.*(\s{1,}|\_{1,})pseudogene$", re.IGNORECASE)
        trnas, genes, pseudogenes, others = list(), list(), list(), list()
        for gt in self._gene_type_read_lengths.keys():
            if re.match(tRNA_regex, gt):
                trnas.append(gt)
            elif re.match(pseudogene_regex, gt):
                pseudogenes.append(gt)
            elif re.match(gene_regex, gt):
                genes.append(gt)
            else:
                others.append(gt)
        gene_types = [
            sorted(others, key=str.lower),
            sorted(trnas, key=str.lower),
            sorted(genes, key=str.lower),
            sorted(pseudogenes, key=str.lower),
        ]
        return list(chain(*gene_types))

    def plot(self) -> None:
        gene_type = "all"
        gt_plot = figure(
            title=f"Length distribution of aligned reads: {gene_type}",
            x_axis_label="Read lengths",
            y_axis_label="Counts",
            width=1200,
            height=900,
        )
        legends, col_pickers, data_source = list(), list(), list()
        colors = sample(Turbo256, len(self.samples))
        for i, ds in enumerate(self.source_data[gene_type]):
            i_cds = ColumnDataSource(data={"rl": ds[0], "counts": ds[1]})
            linep = gt_plot.line(
                x="rl",
                y="counts",
                source=i_cds,
                line_width=2.5,
                line_alpha=0.8,
                color=colors[i],
            )
            data_source.append(i_cds)
            sample_name = re.sub(r"(\_{1,}|\-{1,}|\.{1,})", " ", self.samples[i])
            # color picker
            picker = ColorPicker(title=sample_name, color=colors[i])
            picker.js_link("color", linep.glyph, "line_color")
            col_pickers.append(picker)
            # legend
            legends.append((sample_name, [linep]))
        # legend format
        legend = Legend(items=legends, orientation="horizontal", ncols=self.nsamples)
        legend.click_policy = "hide"
        gt_plot.add_layout(legend, "below")
        # title
        gt_plot.title.align = "center"
        gt_plot.title.text_font_size = "16pt"
        gt_plot.title.text_font_style = "bold"
        # limit zoom limit
        gt_plot.x_range.min_interval = 1
        # axis labels
        # X
        gt_plot.xaxis.major_label_text_font_size = "12pt"
        gt_plot.xaxis.axis_label_text_font_size = "14pt"
        gt_plot.xaxis.axis_label_text_font_style = "bold"
        # Y
        gt_plot.yaxis.major_label_text_font_size = "12pt"
        gt_plot.yaxis.axis_label_text_font_size = "14pt"
        gt_plot.yaxis.axis_label_text_font_style = "bold"
        # add comma to counts on Y axis
        gt_plot.yaxis.formatter = NumeralTickFormatter(format="0,0")
        # gene type selector
        gt_select = Select(
            title="Select gene type",
            value="all",
            options=self._group_n_sort_gene_types(),
        )
        # custom js code
        js_code = CustomJS(
            args=dict(
                source_map=self.source_data,
                data_source=data_source,
                title=gt_plot.title,
            ),
            code="""
            let gt_source_map = source_map[this.value];
            title.text='Length distribution of aligned reads: '+ this.value;
            console.log(title.text);
            for(let i=0; i<Object.entries(data_source).length; i++){
                data_source[i].data['rl'] = gt_source_map[i][0];
                data_source[i].data['counts'] = gt_source_map[i][1];
                data_source[i].change.emit();
            }
        """,
        )
        gt_select.js_on_change("value", js_code)
        if Path(self.output_file).exists():
            logging.warning(f"Over writing {self.output_file}")
        # organize per sample color picker according to self.nsamples
        if self.nsamples > 1:
            col_picks = [
                row(*col_pickers[i : i + self.nsamples])
                for i in range(0, len(col_pickers), self.nsamples)
            ]
            save(
                row(gt_plot, column(gt_select, column(*col_picks))),
                filename=self.output_file,
                title="Gene type read length distribution",
            )
        else:
            save(
                row(gt_plot, column(gt_select, column(*col_pickers))),
                filename=self.output_file,
                title="Gene type read length distribution",
            )
