kraken
------

This group contains a single command that takes NCBI taxonomy .dmp files and `kraken2`_ report files as inputs and generates a combined csv output with taxonomic information and read counts.

collect-reports
^^^^^^^^^^^^^^^

Aggregate `kraken2`_ species classification report. 

Given NCBI taxonomy database `nodes.dmp`_ file, `names.dmp`_ file, either a directory with Kraken2 classification results or a list of Kraken2 report files, aggregate all classification results into one file.

**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value"

    "--nodes", "NCBI taxonomy database `nodes.dmp`_ file", "|tick|"
    "--names", "NCBI taxonomy database `names.dmp`_ file", "|tick|"
    "--out-csv", "Output csv file name", "|tick|"
    "--kraken-dir", "Directory containing `kraken2`_ classification reports. Either --kraken-dir or reports as space separated arguments MUST be provided", "|cross|"
    "--pattern", "Kraken report naming pattern. Used only if `--kraken-dir` is provided", "|cross|", "`*.txt`"

**Usage**

using ``--kraken-dir``

.. code-block:: bash

    # using --kraken-dir
    ngs-statter collect-reports --kraken-dir path/to/kraken_dir --out-csv path/to/reports.csv

using space separated kraken report files as arguments

.. code-block:: bash

    # using space separated kraken report files as arguments
    ngs-statter collect-reports --out-csv path/to/reports.csv path/to/report1.txt path/to/report2.txt 


.. _`kraken2`: https://github.com/DerrickWood/kraken2
.. _`nodes.dmp`: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt
.. _`names.dmp`: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt