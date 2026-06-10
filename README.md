
## NGS statter

A package to compute read/alignment statistics per fastq/bam files. This package is used as a part of an NGS analysis pipeline to generate intermediate QC results and final alignment statistics for all samples.

**Parsers**  
- [x] Fastq parser
- [x] GFF parser
- [x] BAM parser
- [x] [kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) report aggregator
- [x] [fastp](https://github.com/OpenGene/fastp) fastp trimming stats aggregator
- [x] [STAR aligner](https://github.com/alexdobin/STAR) alignment stats to json
- [x] parser fix UMI position in flexbar fastqs

**Plotting functions**  
- [x] Fastq read length distribution plot
- [ ] Overall aligned read length distribution plot
- [x] Gene type specific read length distribution plot
- [ ] Alignment statistics

## Commands
use ``` statter -h ``` for a list of available helper commands.

Commands are grouped based on the input as follows:


#### `alignment` 

Parse alignment (.bam) files and generate base statistics
| command | description | usage |
|:----------:|:----------|:----------|
| bam    |Compute basic alignment statistics for a generic BAM file| `ngs-statter bam -h`|
| STAR|Compute alignment statistics for a STAR aligned BAM file and output as JSON|`ngs-statter STAR -h`|

#### `crosslink`
Parse output crosslink (.bed) files from [Shoji](https://github.com/EMBL-Hentze-group/Shoji) or [htseq-clip](https://github.com/EMBL-Hentze-group/htseq-clip)

| command | description | usage |
|:----------:|:----------|:----------|
| csv-meta-example|Print example CSV metadata file for crosslink plotting to console| `ngs-statter csv-meta-example -h`|
|count-crosslinks | Count crosslinking sites over secondary structure/primary motif regions|`ngs-statter count-crosslinks -h`|
|crosslink-line-plot |Plot crosslinking sites over secondary structure/primary motif regions as line plots|`ngs-statter crosslink-line-plot -h`|
|crosslink-heatmap | Plot crosslinking sites over secondary structure/primary motif regions as heatmaps|`ngs-statter crosslink-heatmap -h`|

#### `fastq`
Parse fastq files to compute and plot read-length distributions

| command | description | usage |
|:----------:|:----------|:----------|
| parse-read-length|Compute fastq read length distribution| `ngs-statter parse-read-length -h`|
| plot-read-length|Read in json formatted output files (from `parse-read-length`) and output html plot| `ngs-statter plot-read-length -h`|

#### `flexbar`
Fix the position of UMIs in flexbar output fastq files
| command | description | usage |
|:----------:|:----------|:----------|
| fix-header |Fix flexbar UMI headers| `ngs-statter fix-header  -h`|

### `genetype`

Given a gff3 formatted gene annotation file and an alignment (.bam) file, compute and plot read-length distribution of reads aligning to gene types (protein coding, lncRNA,...).
| command | description | usage |
|:----------:|:----------|:----------|
| parse-gene-type-read-length |Compute gene type read length distributions| `ngs-statter parse-gene-type-read-length  -h`|
| plot-gene-type-read-length |Plot gene type read length distributions using outputs from `parse-gene-type-read-length`| `ngs-statter plot-gene-type-read-length  -h`|

### [`kraken`](https://github.com/DerrickWood/kraken2)
Given NCBI taxonomy files and a set of Kraken2 species classification reports, merge the reports into one single output.
| command | description | usage |
|:----------:|:----------|:----------|
| collect-reports |Collect kraken2 reports into a single file| `ngs-statter collect-reports -h`|

### sample

| command | description | usage |
|:----------:|:----------|:----------|
|sample-stats|Collect statistics from various steps in a workflow (trimming, alignment,...) and compile into a single json file| `ngs-statter sample-stats -h`|
|compile-stats|Collect stats from multiple samples (output from `sample-stats`) into a single tsv file| `ngs-statter compile-stats -h`|