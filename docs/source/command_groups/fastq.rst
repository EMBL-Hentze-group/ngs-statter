fastq
-----
Commands under this group takes fastq files as input, compute read length distribution and plot them as interactive line plots.

.. _fq-read-len:
parse-read-length
^^^^^^^^^^^^^^^^^

For a given fastq file, compute read lengths and write the distribution to json format

**Options**

.. csv-table::
    :widths: 10 70 10
    :header: "option", "description", "required"

    "--fastq", "Fastq file (supports .gz files)", "|tick|"
    "--json", "File name for json formatted output file", "|tick|"

**Usage**

.. code-block:: bash

    ngs-statter parse-read-length --fastq path/to/fastq.gz --json path/to/output.json

plot-read-length
^^^^^^^^^^^^^^^^
Given a list of json files containing read length distribution (generated using :ref:`fq-read-len`), plot the read length distribution as an interactive (html) line plot.


**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value"

    "--json-dir", "Directory containing json formatted read_length files (see :ref:`fq-read-len`). Either `--json-dir` or json formatted files as space separated arguments MUST be provided", "|cross|"
    "--html", "Output file name", "|tick|"
    "--min-read-length", "Minimum read length that will be kept", "|cross|", 15
    "--pattern", "If `--json-dir` is provided, use this pattern to match files", "|cross|", "`*.json`"
    "--title", "Plot title", "|cross|", "Read length distribution after adapter trimming"

**Usage**

using ``--json-dir``

.. code-block:: bash

    # using --json-dir
    ngs-statter plot-read-length --json-dir path/to/json/dir --html path/to/output.html 

using space separated json files as arguments

.. code-block:: bash

    # using space separated json files as arguments
    ngs-statter plot-read-length --html path/to/output.html path/to/fq1.json path/to/fq2.json 