
## Mad statter

A package to compute read/alignment statistics per fastq/bam file and parsers for niche formats. This will be eventually part of the *seq data analysis workflow. 

**Parsers**  
- [x] Fastq parser
- [x] GFF parser
- [x] BAM parser
- [x] single end unmapped fastq parser
- [ ] paired end unmapped fastq parser
- [x] [ganon](https://pirovc.github.io/ganon/) unclassified single end reads parser
- [x] [kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) report aggregator

**Plotting functions**  
- [x] Fastq read length distribution plot
- [ ] Overall aligned read length distribution plot
- [x] Gene type specific read length distribution plot
- [ ] Alignment statistics
- [x] Motif logos from MMseqs2 clusters
- [x] Ganon contamination and unclassified read length distribution plot

## How to pull the container image from this repo

* Go to [container registry](https://git.embl.de/grp-hentze/workflows/containers/mad_statter/container_registry),  
* Click on container name (mad_statter here), this will show all the tags available in this registry
* Choose the latest published tag

```bash
$ cd /path/to/Singularity_contaier/dir
$ singularity pull --docker-login docker://registry.git.embl.de/grp-hentze/workflows/containers/mad_statter:<container_tag> # This will prompt for git user name and password
Enter Docker Username:<git.embl.de user name>
Enter Docker Password:<password>
# The image will be pulled into this directory as mad_statter_<TAG>.sif
$ ls -alh
```

## Scripts

- read_length_parser: for a given fastq file, dump read length distribution in json
- read_length_plotter: generate interactive html plots from fastq read length distributions
- gene_type_read_length_parser: for a given gff3 file and bam file, dump read length distribution per aligned gene type to json  
- gene_type_read_length_plotter: generate interactive html plots from bam read length distributions
- unmapped_fastq_single_end: split unmapped reads to multimapper unmapped and other unmapped reads
- ganon_unc_to_fastq: given a ganon classification file and corresponding fastq file, create a new fastq file with all ganon unclassified reads
- ganon_read_length_plotter: for ganon classification result plot contamination vs unclassified read length distribution
- ganon_abundance_aggregator: aggregate ganon classification report from multiple samples into a single file

