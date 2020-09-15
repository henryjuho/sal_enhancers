# Construction of a training set for BQSR and VQSR

We used a training set consisting of variants called by both gatk and samtools, and also used on available SNP chips.

GATK calling:

```
mkdir /scratch/project_2002047/barson_reseq/
mkdir -p /scratch/project_2002047/barson_reseq/training/gatk
mkdir -p /scratch/project_2002047/barson_reseq/training/samtools

python haplotype_caller.py -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -bam_dir /scratch/project_2002047/barson_mapping_v2/merged_bams/ -out_dir /scratch/project_2002047/barson_reseq/training/gatk/

mkdir /scratch/project_2002047/barson_reseq/training/gatk_autosomal/
ls /scratch/project_2002047/barson_reseq/training/gatk/*vcf | python reduce_all_gvcf.py -out /scratch/project_2002047/barson_reseq/training/gatk_autosomal/
ls /scratch/project_2002047/barson_reseq/training/gatk_autosomal/*vcf | python index_g_vcf.py
ls /scratch/project_2002047/barson_reseq/training/gatk_autosomal/*vcf | python combine_gvcfs.py -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -out /scratch/project_2002047/barson_reseq/training/gatk_autosomal/salsal31

mkdir /scratch/project_2002047/barson_reseq/training/gatk_calls
ls /scratch/project_2002047/barson_reseq/training/gatk_autosomal/salsal31_NC_0273*vcf | python genotypeGVCFs.py -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -out /scratch/project_2002047/barson_reseq/training/gatk_calls/
ls /scratch/project_2002047/barson_reseq/training/gatk_calls/*vcf | python merge_vcfs.py -out /scratch/project_2002047/barson_reseq/training/gatk_calls/salsal31.gatk.raw.snps.indels.vcf
gatk SelectVariants -R /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -V /scratch/project_2002047/barson_reseq/training/gatk_calls/salsal31.gatk.raw.snps.indels.vcf --select-type-to-include SNP --exclude-non-variants -O /scratch/project_2002047/barson_reseq/training/salsal31.gatk.raw.snps.vcf
```
SAMtools calling

```
ls /scratch/project_2002047/barson_mapping_v2/merged_bams/*bam > bam_list.txt
python samtools_calling.py -bam_list bam_list.txt -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -out /scratch/project_2002047/barson_reseq/training/samtools/salsal31
ls /scratch/project_2002047/barson_reseq/training/samtools/salsal31.samtools.allsites.NC_0273*vcf | python extract_chromo_snps.py -ref /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna
cp /scratch/project_2002047/barson_reseq/training/samtools/salsal31.samtools.raw.snps.vcf* /scratch/project_2002047/barson_reseq/training/
```
Concordance 

```
cd /scratch/project_2002047/barson_reseq/training/
gatk SelectVariants -R /scratch/project_2002047/sal_reseq_v_mapping/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -V salsal31.gatk.raw.snps.vcf --concordance salsal31.samtools.raw.snps.vcf -O salsal31.gatk.samtools.consensus.raw.snps.vcf
```


Getting bed file of SNPs on chips:

```
cat 60k_200k_chip_pos_1based.txt | python chip2bed.py | sort -k1,1 -k2,2n | bgzip -c > chip_snps.bed.gz
tabix -p bed chip_snps.bed.gz 
```
Intersecting called SNPs with chip SNPs.

```
bedtools intersect -header -a salsal31.gatk.samtools.consensus.raw.snps.vcf -b ~/sal_enhancers/training_set/chip_snps.bed.gz | bgzip -c > salsal31.training.snps.vcf.gz
tabix -pvcf salsal31.training.snps.vcf.gz
```


| Stage        | Number of SNPs |
|:-------------|:--------------:|
| gatk         | 10334809       |
| samtools     | 9658669        |
| intersect    | 7844181        |
| chip         | 216324         |