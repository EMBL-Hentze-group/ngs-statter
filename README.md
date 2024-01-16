
## Mad statter

A package to compute read/alignment statistics per fastq/bam file and parsers for niche formats. This will be eventually part of the *seq data analysis workflow. 

**Parsers**  
- [x] Fastq parser
- [x] GFF parser
- [x] BAM parser
- [x] single end unmapped fastq parser
- [ ] paired end unmapped fastq parser
- [x] [ganon](https://pirovc.github.io/ganon/) unclassified single end reads parser

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

