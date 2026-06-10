sample
------

Commands under this group are meant to be used as a part of an NGS data processing pipeline, to generate read statistics at different stages (trimming, pre-processing, alignment etc) of the pipeline, and finally combine them into a single tsv file for all the samples in the pipeline, for the final report.

.. _sample-stats-ref:
sample-stats
^^^^^^^^^^^^
Aggregate read statistics from various steps in the workflow for a single sample into a json file.

**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value"

    "--sample-name", "Sample name", "|tick|"
    "--raw-reads", "Path to `seqkit stats`_ output for raw reads", "|tick|"
    "--first-trim", "Path to `seqkit stats`_ output for reads after first trimming step", "|tick|"
    "--second-trim", "Path to `seqkit stats`_ output for reads after second trimming step", "|cross|"
    "--rRNA-mapped", "Path to `seqkit stats`_ output for reads mapped to rRNA reference", "|cross|"
    "--rRNA-free", "Path to `seqkit stats`_ output for reads not mapped to rRNA reference", "|cross|"
    "--align-stats", "Path to bam statistics json file generated using :ref:`STAR-ref` command", "|tick|"
    "--dedup-stats", "Path to bam statistics json file after UMI deduplication generated using :ref:`STAR-ref` command", "|cross|"
    "--kraken2-report", "Path to `kraken2`_ report file", "|cross|"
    "--out-json", "Output json file for collected statistics", "|tick|"

**Usage**

.. code-block:: bash

    ngs-statter sample-stats --sample-name SAMPLE_NAME --raw-reads path/to/raw_reads_stats.txt --first-trim path/to/first_trim_stats.txt --align-stats path/to/align_stats.json --out-json path/to/sample_stats.json

compile-stats
^^^^^^^^^^^^
Aggregate read statistics for multiple samples (generated using :ref:`sample-stats-ref` command) into a single tsv file for the final report. 

**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value"

    "--output", "Output file name for all sample statistics", "|tick|"

**Usage**

.. code-block:: bash

    ngs-statter compile-stats --output path/to/compiled_stats.tsv path/to/sample1_stats.json path/to/sample2_stats.json

.. _`kraken2`: https://github.com/DerrickWood/kraken2
.. _`seqkit stats`: https://bioinf.shenwei.me/seqkit/usage/#stats