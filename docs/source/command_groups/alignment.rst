alignment
---------

The following tools are grouped under alignment as they take alignment files (.bam) as input.

bam
^^^

Given a bam file as input, compute the following general alignment statistics and output them as json. The output will contain the following fields:

.. csv-table::
    :header: "field", "description"

    "Input reads","total number of reads in the BAM file"                                                
    "Mapped","total number of reads mapped to the reference genome"                                    
    "Unmapped","total number of reads that are not mapped"

**Options**


.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value" 

    "--bam", "Path to the input bam file (MUST be co-ordinate sorted and indexed)", "|tick|"
    "--out-json", "Output file to write json formatted data", "|tick|"
    "--min-q","Minimum alignment quality", "|cross|", "0"

**Usage**

.. code-block:: bash

    ngs-statter bam --bam path/to/alignment.bam --out-json path/to/output.json --min-q 0
    
.. _STAR-ref:
STAR
^^^^

Compute alignment statistics for an alignment file generated using `STAR aligner`_ and output as json. The output will contain the following fields:

.. csv-table::
    :header: "field", "description"

    "Reads for mapping","total number of reads in the BAM file"                                                
    "Mapped: Total","total number of reads mapped to the reference genome"
    "Mapped: Uniquely mapped reads", "total number of reads mapped to a unique location in the reference genome"
    "Mapped: Multimapped reads", "total number of reads mapped to multiple locations in the reference genome"
    "Mapped: PCR duplicate reads", "total number of mapped reads marked as PCR duplicates in the BAM file, if the BAM file is marked for duplicates, see `samtools markdup`_"
    "Mapped: Unique reads", "total number of reads mapped to a unique location in the reference genome and not marked as PCR duplicates in the BAM file, if the BAM file is marked for duplicates, see `samtools markdup`_"
    "Unmapped: Total","total number of reads that are not mapped to the reference genome"
    "Unmapped: mapped to too many loci", "total number of reads that marked as unmapped as they are mapped to too many locations in the reference genome"
    "Unmapped: no seed/windows", "total number of reads that are marked as unmapped as they do not have a seed region that can be mapped to the reference genome"
    "Unmapped: too many mismatches", "total number of reads that are marked as unmapped as they have too many mismatches compared to the reference genome"
    "Unmapped: too short", "total number of reads that are marked as unmapped as the seed regions are too short to be mapped to the reference genome"
    "Unmapped: paired-end mate","for paired end reads, total number of reads that are marked as unmapped as their paired end mate is mapped to the reference genome"

**Options**


.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value" 

    "--bam", "Path to the input bam file (MUST be co-ordinate sorted and indexed)", "|tick|"
    "--out-json", "Output file to write json formatted data", "|tick|"
    "--min-q","Minimum alignment quality", "|cross|", "0"

**Usage**

.. code-block:: bash

    ngs-statter STAR --bam path/to/star_alignment.bam --out-json path/to/output.json --min-q 0
    

.. _`STAR aligner`: https://github.com/alexdobin/STAR
.. _`samtools markdup`: https://www.htslib.org/doc/samtools-markdup.html