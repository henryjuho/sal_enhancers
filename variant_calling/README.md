# SNP calling for Barson et al. 2015 re-sequencing data

## BQSR

SNP calling was performed using GATK (version 4.1.3.0). First off bam files underwent BQSR (base quality score recalibration), using the sites called by both GATK and SAMtools and found on the salmon SNP chip as training data (see: ) :

```
mkdir /scratch/project_2002047/barson_mapping_v2/bqsr_bams
ls /scratch/project_2002047/barson_mapping_v2/merged_bams/*bam | python BQSR.py -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -train_vcf /scratch/project_2002047/barson_reseq/training/salsal31.training.snps.vcf.gz -out_dir /scratch/project_2002047/barson_mapping_v2/bqsr_bams/
```

## SNP calling

GATKs haplotype caller was then used to create ```g.vcf``` files per sample per chromosome. Scaffolds were exlcuded. This uses a modified wrapper script to the one used in creating the training set. There seems to be a bug in combineGVCFs, where if scaffolds are listed in the header it freezes even though there is only one chromosome in each ```.g.vcf``` file, so trimmed the headers.

```
mkdir /scratch/project_2002047/barson_reseq/main_calling
ls /scratch/project_2002047/barson_mapping_v2/bqsr_bams/*bam | python haplotype_caller_II.py -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -out_dir /scratch/project_2002047/barson_reseq/main_calling/

mkdir /scratch/project_2002047/barson_reseq/trimmed_head_gvcf
python trim_gvcf_contig_header.py -in_dir /scratch/project_2002047/barson_reseq/main_calling/ -out_dir /scratch/project_2002047/barson_reseq/trimmed_head_gvcf/
```

Individual ```.g.vcf``` files were then merged into mutlisample ones:

```
mkdir /scratch/project_2002047/barson_reseq/chromosomal_allsites
python combine_gvcfs_II.py -in_dir /scratch/project_2002047/barson_reseq/trimmed_head_gvcf/ -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -out_dir /scratch/project_2002047/barson_reseq/chromosomal_allsites/
```

Variants were then called, and autosomal ```.vcf``` files merged.

```
mkdir /scratch/project_2002047/barson_reseq/genotyped_variants
ls /scratch/project_2002047/barson_reseq/chromosomal_allsites/salsal_31.NC_0273*vcf | python ../training_set/genotypeGVCFs.py -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -out /scratch/project_2002047/barson_reseq/genotyped_variants/
ls /scratch/project_2002047/barson_reseq/genotyped_variants/*vcf | python ../training_set/merge_vcfs.py -out /scratch/project_2002047/barson_reseq/genotyped_variants/salsal_31.autosomes.raw.snps.indels.vcf
```

## VQSR and filtering

We performed VQSR (variant quality score recalibration) using the same training set as BQSR.

```
mkdir /scratch/project_2002047/barson_reseq/vqsr
ls /scratch/project_2002047/barson_reseq/genotyped_variants/*vcf | python ../training_set/merge_vcfs.py -out /scratch/project_2002047/barson_reseq/genotyped_variants/salsal_31.autosomes.raw.snps.indels.vcf
```

| Tranche: | 100.0  | 99.9   | 99.5   | 99.0   | 98.0   | 95.0   | 90.0  |
|:--------:|:------:|:------:|:------:|:------:|:------:|:------:|:-----:|
| N SNPs:  |9736064 |6797538 | **5701538** |5305874 |5009797 |4862164 |3754860|

We retained variants passing the 99.5 tranche cut off and then applied repeat filters, depth filters (twice and half the mean depth of 8X), and removed SNPs with missing genotypes.

```
mkdir /scratch/project_2002047/barson_reseq/post_vqsr
gatk --java-options "-Xmx14G" VariantFiltration -R /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -V /scratch/project_2002047/barson_reseq/vqsr/salsal_31.autosomes.raw.snps.indels.recalibrated.filtered_t99.5.pass.vcf -O /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.raw.snps.indels.recalibrated.filtered_t99.5.pass.rmarked.vcf --mask /scratch/tuyida/bartonhe/sal_ref/salmo_salar_repeats.bed.gz --mask-name REPEAT
 gatk --java-options "-Xmx14G" SelectVariants -R /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -V /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.raw.snps.indels.recalibrated.filtered_t99.5.pass.rmarked.vcf -O /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.raw.snps.indels.recalibrated.filtered_t99.5.pass.rfiltered.biallelic.vcf --exclude-filtered -restrict-alleles-to BIALLELIC

python ~/sal_bal_sel/training_set/depth_filter.py -vcf /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.raw.snps.indels.recalibrated.filtered_t99.5.pass.rfiltered.biallelic.vcf -mean_depth 8 -N 31
Depth range: 4.0 - 16.0
4837611 variants passed, 241336 variants failed
Output written to /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.raw.snps.indels.recalibrated.filtered_t99.5.pass.rfiltered.biallelic.dpfiltered.vcf

python filter_uncalled_snps.py -vcf /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.raw.snps.indels.recalibrated.filtered_t99.5.pass.rfiltered.biallelic.dpfiltered.vcf | bgzip -c > /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.t99_5_snps.allfilters.vcf.gz
tabix -pvcf  /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.t99_5_snps.allfilters.vcf.gz
```

Summary of filtering:


| Stage | Number of SNPs |
|:------|:--------------:|
| gatk raw  | 9736064  |
| VQSR (99.5) |5701538 |
| repeats and biallelic |5078947 |
| depth | 4837611  |
| missing calls | 3723849 |

## Callable sites

A fasta file of callable site in the genome was created (i.e. those monomorphic sites that passed all the filters).

```
mkdir /scratch/project_2002047/barson_reseq/callable_sites
ls /scratch/project_2002047/barson_reseq/chromosomal_allsites/*vcf | while read i; do bgzip $i; done
ls /scratch/project_2002047/barson_reseq/chromosomal_allsites/*vcf.gz | while read i; do tabix -pvcf $i; done

ls /scratch/project_2002047/barson_reseq/chromosomal_allsites/*vcf.gz |  python callable_sites_parallel.py -bed /scratch/tuyida/bartonhe/sal_ref/salmo_salar_repeats.bed.gz -chr_bed /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_autosomes.bed -pol /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -out_pre /scratch/project_2002047/barson_reseq/callable_sites/salsal_31.callable
cat /scratch/project_2002047/barson_reseq/callable_sites/salsal_31.callable_falist.txt | python check_fasta_chunks.py

cat /scratch/project_2002047/barson_reseq/callable_sites/salsal_31.callable_falist.txt | python fa_cat.py > /scratch/project_2002047/barson_reseq/callable_sites/salsal_31.callable.fa
samtools faidx /scratch/project_2002047/barson_reseq/callable_sites/salsal_31.callable.fa

python check_fasta.py /scratch/project_2002047/barson_reseq/callable_sites/salsal_31.callable.fa /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna
```