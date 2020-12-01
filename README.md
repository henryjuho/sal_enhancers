# DFE analysis pipeline from Verta et al. (2020 - preprint)

This repository contains the code and command lines for the re-sequencing calling pipeline and DFE analysis presented in 
Verta et al. (2020 - preprint) (<https://www.biorxiv.org/content/10.1101/2020.11.20.389684v1.abstract>). The different 
analysis and preparation steps are contained in their own subdirectories which are outlined below in rough sequential 
order.

## Contents 

[annotation/](annotation/): containes the code for downloading the reference genome and coordinate files from NCBI, as well as code to calculate coordinates for sites without coordinates in the GFF file.

[genome_alignment/](genome_alignment/): pipeline for 9 species whole genome alignment - used for polarising SNPs and calculating divergence in order to calculate alpha.

[read_mapping/](read_mapping/): steps for downloading the raw reads from Barson et al. (2015) and mapping to the reference genome.

[training_set/](training_set/): initial SAMtools and GATK calling to produce a training set for the main GATK calling.

[variant_calling/](variant_calling/): GATK SNP calling pipeline and generation of callable sites FASTA file.

[sfs/](sfs/): Pipeline to extract frequency data and number of callable sites for different sequence contexts per gene.

[summary_stats/](summary_stats/): Calculating nucleotide diversity, theta and Tajima's D from SFS data.

[dfe/](dfe/): Estimating the distribution of fitness effects using anavar.

[divergence/](divergence/): Estimating divergence from the genome alignment using APE and estimating alpha with the DFE and divergence estimates.

## A note on the cluster

A lot of the pipeline scripts write and submit jobs on SLURM based CSC computer cluster Puhti (<https://docs.csc.fi/computing/overview/>).
This behaviour can be changed by substituting the line ```from qsub import q_sub``` with ```from qsub import q_write as q_sub``` at the top
of the offending python scripts.

## Python requirements

The scripts also make use of a number of python modules:

```anavar_utils```: <https://github.com/henryjuho/anavar_utils> - some code to make working with anavar control files easier.

```python_qsub_wrapper```: <https://github.com/henryjuho/python_qsub_wrapper> - some code to write and submit batch jobs from within the python scripts.

```sfs_utils```: <https://github.com/henryjuho/sfs_utils> - code for extracting frequency data from VCF files

```pysam```: <https://pysam.readthedocs.io/en/latest/api.html> - module for working with VCF, BED, FASTA and other files in python.

## Other directories

There are also the following directories that contain pipelines that did not end up in the paper:

[homeoblock_alignments/](homeoblock_alignments/): Pipeline for align homeoblocks in the salmon, not used for any analysis.

[coverage_impact/](coverage_impact/): A sanity check to see if including low coverage individuals to maximise sample size and number of SNPs was skewing the estimated DFE. 
