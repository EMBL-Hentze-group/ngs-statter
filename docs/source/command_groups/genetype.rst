genetype
--------

The commands under this group takes a gene annotation file in gff format and an alignment file in bam format as inputs, computes and plots the distribution of read lengths over different gene types (protein coding, lncRNA, snoRNA etc) as interactive line plots.

.. _gene-type-length-ref:
parse-gene-type-read-length
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Given a gene annotation file in gff format and an alignment file in bam format, compute the distribution of read lengths over different gene types (protein coding, lncRNA, snoRNA etc), and write the output in json format.

**Options**

.. csv-table::
    :widths: 10 70 5 15
    :header: "option", "description", "required", "default value" 

    "--gff3", "Gene annotation file in gff format (supports .gz files)", "|tick|"
    "--bam", "Alignment file in bam format (MUST be co-ordinate sorted and indexed)", "|tick|"
    "--out-json", "Output file to write json formatted data", "|tick|"
    "--gene-features","Features to parse from GFF3 file, MUST be supplied as: ``--gene-features gene --gene-features tRNA...``", "|cross|", "``--gene-features gene --gene-features tRNA``"
    "--gene-id", "Gene id attribute in GFF3 attribute column", "|cross|", "gene_id"
    "--gene-type", "Gene type attribute in GFF3 attribute column", "|cross|", "gene_type"
    "--min-q","Minimum alignment quality", "|cross|", "0"
    "--ignore-duplicate", "Flag to ignore PCR duplicate reads (after `samtools markdup`_)"

**Usage**

.. code-block:: bash

    ngs-statter parse-gene-type-read-length --gff3 path/to/annotation.gff3.gz --bam path/to/alignment.bam --out-json path/to/output.json

plot-gene-type-read-length
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Given a list of json files containing read length distribution over gene types (generated using :ref:`gene-type-length-ref`), plot the read length distribution as an interactive (html) line plot.

**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value"

    "--json-dir", "Directory containing json formatted read length per gene type files (see :ref:`gene-type-length-ref`). Either `--json-dir` or json formatted files as space separated arguments MUST be provided", "|cross|"
    "--html", "Output file name", "|tick|"
    "--pattern", "If `--json-dir` is provided, use this pattern to match files", "|cross|", "`*.json`"

**Usage**

using ``--json-dir``

.. code-block:: bash

    # using --json-dir
    ngs-statter plot-gene-type-read-length --json-dir path/to/json_dir --html path/to/output.html 

using space separated json files as arguments

.. code-block:: bash

    # using space separated json files as arguments
    ngs-statter plot-gene-type-read-length --html path/to/output.html path/to/bam1.json path/to/bam2.json 

.. _`samtools markdup`: https://www.htslib.org/doc/samtools-markdup.html