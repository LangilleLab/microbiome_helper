Microbiome Helper
=================

An assortment of scripts to help process and automate various microbiome and metagenomic bioinformatic tools.

Including workflows (see below) for analyzing 16S and shotgun metagenomic data. These workflows can be run on an [Ubuntu virtual box image](https://github.com/mlangill/microbiome_helper/wiki/MicrobiomeHelper-Virtual-Box-image) for those who don't want to set up their own environment.

## Table of contents:

* [Requirements](https://github.com/mlangill/microbiome_helper/wiki/Requirements)
* [Brief description of scripts](https://github.com/mlangill/microbiome_helper/wiki/Brief-description-of-scripts)
* [Virtual Box image](https://github.com/mlangill/microbiome_helper/wiki/MicrobiomeHelper-Virtual-Box-image)

### Metagenomic resources
   * [Metagenomic standard operating procedure](https://github.com/mlangill/microbiome_helper/wiki/Metagenomic-standard-operating-procedure)

**More detail on particular steps:**
   * [Stitch reads](https://github.com/mlangill/microbiome_helper/wiki/Stitch-reads)
   * [Sequence QC](https://github.com/mlangill/microbiome_helper/wiki/Sequence-QC)
   * [Screen out human sequences](https://github.com/mlangill/microbiome_helper/wiki/Screen-out-human-sequences)
   * [Taxonomic composition](https://github.com/mlangill/microbiome_helper/wiki/Taxonomic-composition)
   * [Functional profiling](https://github.com/mlangill/microbiome_helper/wiki/Functional-profiling)

### 16S resources
   * [16S standard operating procedure](https://github.com/mlangill/microbiome_helper/wiki/16S-standard-operating-procedure)
   * [16S tutorial](https://github.com/mlangill/microbiome_helper/wiki/16S-tutorial)

**More detail on particular steps:**
   * [Stitch reads](https://github.com/mlangill/microbiome_helper/wiki/Stitch-reads)
   * [Sequence QC](https://github.com/mlangill/microbiome_helper/wiki/Sequence-QC)
   * [Remove low quality reads](https://github.com/mlangill/microbiome_helper/wiki/Remove-low-quality-reads)
   * [Remove chimeric reads](https://github.com/mlangill/microbiome_helper/wiki/Remove-chimeric-reads)
   * [OTU picking](https://github.com/mlangill/microbiome_helper/wiki/OTU-picking)
   * [Alpha and beta diversity](https://github.com/mlangill/microbiome_helper/wiki/Alpha-and-beta-diversity)
   * [Changes to 16S workflow for 18S data](https://github.com/mlangill/microbiome_helper/wiki/Changes-to-16S-workflow-for-18S-data)
   * [PICRUSt workflow](https://github.com/mlangill/microbiome_helper/wiki/PICRUSt-workflow)
   * [Additional QIIME analysis](https://github.com/mlangill/microbiome_helper/wiki/Additional-QIIME-analysis)

### Before getting started:

* These workflows start with raw paired-end MiSeq data in demultiplexed fastq format.

* Most of the commands allow parallel threads to be used with the -p option. Adjust according to your computer hardware abilities.

* These workflows can be completed in roughly 24 hours on a single quad-core desktop when starting with 25 million paired end reads.   


Contact
-------

Questions or comments about this repository can be sent to morgan.g.i.langille@dal.ca


