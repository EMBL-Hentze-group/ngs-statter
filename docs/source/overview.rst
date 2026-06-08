overview
==============

This section provides an overview of the tools available in NGS statter package. The commands in this package are divided into the following groups based on their purpose and input data types.

NGS statter has the following command groups:

alignment
^^^^^^^^^
The following tools are grouped under alignment as they take alignment files (.bam) as input.

crosslink
^^^^^^^^^
Commands under this group takes crosslinking data (.bed) in bed format, either from Shoji or htseq-clip or other similar tools and a bed formatted Region Of Interest (ROI) file as inputs, calculate the crosslink profile over the regions and flank and plot them either as line plots or heatmaps.

fastq
^^^^^
Commands under this group takes fastq files as input, compute read length distribution and plot them as interactive line plots.

flexbar
^^^^^^^
This group contains a single command that takes fastq files from flexbar as input and fix the location of UMIs in the headers.

genetype
^^^^^^^^
The commands under this group takes a gene annotation file in gff format and an alignment file in bam format as inputs, computes and plots the distribution of read lengths over different gene types (protein coding, lncRNA, snoRNA etc) as interactive line plots.

kraken
^^^^^^

This group contains a single command that takes NCBI taxonomy .dmp files and `kraken2`_ report files as inputs and generates a combined csv output with taxonomic information and read counts.

sample
^^^^^^
Commands under this group are meant to be used as a part of an NGS data processing pipeline, to generate read statistics at different stages (trimming, pre-processing, alignment etc) of the pipeline, and finally combine them into a single tsv file for all the samples in the pipeline, for the final report.


.. _`kraken2`: https://github.com/DerrickWood/kraken2
