VCF
---------

Commands under this group takes Crosslink Induced Mutation Sites (CIMS) in VCF format and a bed formatted Region Of Interest (ROI) file as inputs, calculate the alternate allele depth coverage over the regions and flank and plot them as line plots.

CIMS are genereated for R2-CLIP data as a part of the standard analysis pipeline.

.. _vcf-csv-ref:
vcf-meta-example
^^^^^^^^^^^^^^^^
Print example CSV metadata file for CIMS plotting                                                                                                       

**Usage**

.. code-block:: bash

    ngs-statter vcf-meta-example

A format example and additional descriptions will be printed to the console.

parse-alt-allele-depth
^^^^^^^^^^^^^^^^

Given a region of interest (ROI) file in bed format and a metadata CSV file (see :ref:`vcf-csv-ref`), this command counts the coverage of alternate alleles in  the regions specified in the ROI file.

**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value" 

    "--metadata", "CSV metadata file specifying VCF files and sample information (see :ref:`vcf-csv-ref`)", "|tick|"
    "--bed", "BED file specifying secondary structure/ primary motif regions of interests (supports .gz files)", "|tick|"
    "--out-table","Output file to write the aggregated table (always `.parquet`_ format)", "|tick|"
    "--l","5' extension length for regions in BED file", "|cross|", "100"
    "--r","3' extension length for regions in BED file", "|cross|", "100"
    "--unstranded", "If this flag is set, ignore strand information in the BED file and treat all regions as unstranded", "|cross|"
    "--most-5prime", "If bed regions overlap, only keep the most 5' region out of the overlapping regions", "|cross|"
    "--sw", "When plotting smooth crosslink sites using moving average. Use these many adjacent bases to compute moving average", "|cross|", 5
    "--variant-type", "Type of variant to consider: ['all', 'snp', 'mnp', 'indel', 'insertion', 'deletion']", "|cross|", "all"
    "--read-depth", "Minimum read depth to consider a variant", "|cross|", 0
    "--genome-fasta", "Reference genome in FASTA format (for normalization)", "|cross|"
    "--tmpdir", "Temporary directory to use (default: system temp folder)", "|cross|"
    "--threads ", "Number of threads to use", "|cross|", 4

**Usage**

.. code-block:: bash

    ngs-statter parse-alt-allele-depth --metadata path/to/metadata.csv --bed path/to/roi.bed --out-table path/to/output.parquet --l 100 --r 100 --variant-type snp --read-depth 5

plot-alt-allele-depth
^^^^^^^^^^^^^^^^

Given a region of interest (ROI) file in bed format and a metadata CSV file (see :ref:`vcf-csv-ref`), this command counts the coverage of alternate alleles in  the regions specified in the ROI file, and plots alternate allele depth profiles as line plots.

**Options**

.. csv-table::
    :widths: 10 70 10 10
    :header: "option", "description", "required", "default value" 

    "--metadata", "CSV metadata file specifying VCF files and sample information (see :ref:`vcf-csv-ref`)", "|tick|"
    "--bed", "BED file specifying secondary structure/ primary motif regions of interests (supports .gz files)", "|tick|"
    "--out-table","Output file to write the aggregated table (always `.parquet`_ format)", "|tick|"
    "--l","5' extension length for regions in BED file", "|cross|", "100"
    "--r","3' extension length for regions in BED file", "|cross|", "100"
    "--unstranded", "If this flag is set, ignore strand information in the BED file and treat all regions as unstranded", "|cross|"
    "--most-5prime", "If bed regions overlap, only keep the most 5' region out of the overlapping regions", "|cross|"
    "--sw", "When plotting smooth crosslink sites using moving average. Use these many adjacent bases to compute moving average", "|cross|", 5
    "--variant-type", "Type of variant to consider: ['all', 'snp', 'mnp', 'indel', 'insertion', 'deletion']", "|cross|", "all"
    "--read-depth", "Minimum read depth to consider a variant", "|cross|", 0
    "--genome-fasta", "Reference genome in FASTA format (for normalization)", "|cross|"
    "--out-fig", "Output file to write the plot (svg/pdf/png)", "|tick|"
    "--fig-width", "Figure width in centimeters", "|cross|", 30
    "--fig-height", "Figure height in centimeters", "|cross|", 27
    "--xlabel", "X axis label for the plot", "|cross|", "X axis label"
    "--ylabel", "Y axis label for the plot", "|cross|", "Y axis label"
    "--title", "Title for the plot", "|cross|", "Crosslink Induced Mutations: Alternate allele depth"
    "--ymax", "Maximum value for allele depth on y axis (determined from data if not set)", "|cross|"
    "--show-group-mean", "If this flag is set, show the mean allele depth for each group", "|cross|"
    "--errorbar", "Error bar to show for allele depth, see `seaborn errorbar`_ tutorial", "|cross|"
    "--tmpdir", "Temporary directory to use (default: system temp folder)", "|cross|"
    "--threads ", "Number of threads to use", "|cross|", 4

**Usage**

.. code-block:: bash

    ngs-statter plot-alt-allele-depth --metadata path/to/metadata.csv --bed path/to/roi.bed --out-table path/to/output.parquet --out-fig path/to/output.svg --l 100 --r 100 --variant-type snp --read-depth 5


.. _`.parquet`: https://parquet.apache.org/
.. _`seaborn errorbar`: https://seaborn.pydata.org/tutorial/error_bars.html
