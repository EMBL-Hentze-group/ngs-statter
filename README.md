
## Mad statter

A package to compute read/alignment statistics per fastq/bam file and parsers for niche formats. This will be eventually part of the *seq data analysis workflow. 

**Parsers**  
- [x] Fastq parser
- [x] GFF parser
- [x] BAM parser
- [x] [kraken2](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) report aggregator
- [x] [sourmash tax metagenome][https://sourmash.readthedocs.io/en/latest/classifying-signatures.html) Kraken style report aggregator

- [x] [fastp](https://github.com/OpenGene/fastp) fastp trimming stats aggregator
- [x] [STAR aligner](https://github.com/alexdobin/STAR) alignment stats to json
- [x] parser fix UMI position in flexbar fastqs

**Plotting functions**  
- [x] Fastq read length distribution plot
- [ ] Overall aligned read length distribution plot
- [x] Gene type specific read length distribution plot
- [ ] Alignment statistics

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
use ``` statter -h ``` for a list of available helper commands.


For compatibility with existing workflows, the following stand alone scripts are also available:

:warning: This is an incomplete list, kept as a reference. For a complete list of helper commands use ```statter -h```
