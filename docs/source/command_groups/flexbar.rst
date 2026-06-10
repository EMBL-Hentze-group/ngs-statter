flexbar
-------

This group contains a single command that takes fastq files from flexbar as input and fix the location of UMIs in the headers.

fix-header
^^^^^^^^^

flexbar outputs fastq files with UMIs in the header, but the location of UMIs is not optimal. By default UMIs are appended to the comment section of the fastq header like this::

    @VH01545:237:AAHGKVNM5:1:1101:20466:1189 1:N:0:1_GGTTTA
    TTCTAACTAAGGTAAGAGCGGTTCAGCATCCACGCCGCTCTTGCTCTCTGATTTTAGTTGGTATTCTAAGTAA
    +
    IIIIII-IIIIIIIIIII9IIIIIII-II--I-I9II-II------I-999999999II-I9-99-999--99

Note that the UMI (``GGTTTA``) is located at the end of the header, after  ``1:N:0:1``, separated by ``_``.

This command moves the UMI to the end of te sequence identifier, the part before the first space in the header, like this::

    @VH01545:237:AAHGKVNM5:1:1101:20466:1189_GGTTTA 1:N:0:1
    TTCTAACTAAGGTAAGAGCGGTTCAGCATCCACGCCGCTCTTGCTCTCTGATTTTAGTTGGTATTCTAAGTAA
    +
    IIIIII-IIIIIIIIIII9IIIIIII-II--I-I9II-II------I-999999999II-I9-99-999--99

So that tools such as `umi_tools`_ can easily extract the UMIs from the header.

**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value" 

    "--in-fq", "Input FASTQ file with UMI in header comment section, supports .gz files", "|tick|"
    "--out-fq", "Output FASTQ file with fixed UMI position in header, supports .gz files", "|tick|"
    "--separator", "Separator used to separate UMI from the rest of the header", "|cross|", "``_``"

**Usage**

.. code-block:: bash

    ngs-statter fix-header --in-fq path/to/input.fq.gz --out-fq path/to/output.fq.gz

.. _`umi_tools`: https://umi-tools.readthedocs.io/en/latest/QUICK_START.html