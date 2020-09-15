# Downloading and mapping reads

The raw reads from Barson et al 2015 were downloaded from the ENA <https://www.ebi.ac.uk/ena/browser/home>. Download links can be found in [PRJEB10744.txt](PRJEB10744.txt).

```shell script
mkdir /scratch/project_2002047/barson_mapping_v2
mkdir /scratch/project_2002047/barson_mapping_v2/reads

python download_reads.py
```

Reads were cleaned/trimmed using trim galore with cutadapt and then mapped to the salmon reference genome with BWA. Duplicates were marked and read group info added with picard.
The above was performed per fastq pair, the resulting bam files were merged to produce individual bams.

```shell script
mkdir /scratch/project_2002047/barson_mapping_v2/clean_reads
python clean_reads.py -fastq_dir /scratch/project_2002047/barson_mapping_v2/reads/ -out_dir /scratch/project_2002047/barson_mapping_v2/clean_reads/

ls /scratch/project_2002047/barson_mapping_v2/clean_reads/ERR1013*1.fq.gz | python get_read_group_info.py
mkdir /scratch/project_2002047/barson_mapping_v2/bams
python map_reads.py -ref /scratch/project_2002047/sal_reseq/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna -bam_dir /scratch/project_2002047/barson_mapping_v2/bams/ -read_info_file fastq_and_rg.txt

# some failed
ls /scratch/project_2002047/barson_mapping_v2/clean_reads/*1.fq.gz | cut -d '/' -f 6 | cut -d '_' -f 1 > fastqs_bams.txt
ls /scratch/project_2002047/barson_mapping_v2/bams/*.dedup.bam | cut -d '/' -f 6 | cut -d '.' -f 1 >> fastqs_bams.txt
sort fastqs_bams.txt | uniq -c | grep -w '1' | tr -s ' ' | cut -d ' ' -f 3 | while read i; do sbatch --account=project_2002047 /scratch/project_2002047/barson_mapping_v2/bams/${i}_mapping_job.sh ; done

mkdir /scratch/project_2002047/barson_mapping_v2/merged_bams/
python merge_bams.py -bam_in /scratch/project_2002047/barson_mapping_v2/bams/ -bam_out /scratch/project_2002047/barson_mapping_v2/merged_bams/ -info_file fastq_and_rg.txt -ref /scratch/project_2002047/sal_reseq/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna

# rerunning failed merges
ls /scratch/project_2002047/barson_mapping_v2/merged_bams/*bam | cut -d '/' -f 6 | cut -d '.' -f 1 > bams_done.txt
ls /scratch/project_2002047/barson_mapping_v2/merged_bams/*sh | cut -d '/' -f 6 | cut -d '_' -f 1,2,3 >> bams_done.txt
sort bams_done.txt | uniq -c | grep -w '1' | tr -s ' ' | cut -d ' ' -f 3 | while read i; do sbatch --account=project_2002047 /scratch/project_2002047/barson_mapping_v2/merged_bams/${i}_merge_job.sh ; done

ls /scratch/project_2002047/barson_mapping_v2/merged_bams/*.txt | python coverage_summary.py
```
Read coverage in resulting mapped data:

| run accession | sample | mean coverage | median coverage |
|:---|:--:|:--:|:--:|
| ERR1013407 | Alta_12_0001 | 7.520952 | 8 |
| ERR1013495 | Alta_12_0124 | 3.349194 | 3 |
| ERR1013583 | Alta_12_0228 | 7.705269 | 8 |
| ERR1013463 | Arga_12_0075 | 7.728542 | 8 |
| ERR1013471 | Arga_12_0082 | 9.305405 | 10 |
| ERR1013479 | Arga_12_0089 | 13.245793 | 14 |
| ERR1013647 | Jols_13_0001 | 10.631584 | 11 |
| ERR1013639 | Jols_13_0009 | 8.723167 | 9 |
| ERR1013655 | Jols_13_0027 | 7.64258 | 8 |
| ERR1013447 | Nams_12_0024 | 10.828244 | 11 |
| ERR1013455 | Nams_12_0071 | 10.654392 | 11 |
| ERR1013439 | Nams_12_0201 | 6.65685 | 7 |
| ERR1013423 | Naus_12_0014 | 9.373985 | 10 |
| ERR1013415 | Naus_12_0037 | 7.585018 | 8 |
| ERR1013431 | Naus_12_0059 | 12.050208 | 13 |
| ERR1013615 | Repp_12_0007 | 8.360628 | 9 |
| ERR1013623 | Repp_12_0023 | 6.412691 | 6 |
| ERR1013631 | Repp_12_0034 | 12.499478 | 13 |
| ERR1013487 | Uts_11_15 | 7.155256 | 7 |
| ERR1013511 | Uts_11_17 | 19.169955 | 21 |
| ERR1013519 | Uts_11_24 | 9.167601 | 10 |
| ERR1013527 | Uts_11_26 | 5.419665 | 5 |
| ERR1013535 | Uts_11_27 | 6.556109 | 7 |
| ERR1013543 | Uts_11_28 | 6.834338 | 7 |
| ERR1013551 | Uts_11_29 | 4.531277 | 4 |
| ERR1013559 | Uts_11_30 | 6.543761 | 7 |
| ERR1013567 | Uts_11_31 | 6.588743 | 7 |
| ERR1013575 | Uts_11_39 | 2.994131 | 3 |
| ERR1013591 | Uts_11_46 | 7.094451 | 7 |
| ERR1013599 | Uts_11_52 | 3.360987 | 3 |
| ERR1013607 | Uts_11_53 | 8.927331 | 9 |