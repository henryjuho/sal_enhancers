# Alignments of salmon homeoblocks against northern pike

## Genomes

Northern pike and atlantic salmon genomes were downloaded from NCBI. Atlantic salmon was then soft masked (pike already masked).

```shell script
mkdir -p /scratch/project_2002047/sal_enhance/alignments

cd /scratch/project_2002047/sal_enhance/alignments
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/634/155/GCF_004634155.1_Eluc_v4/GCF_004634155.1_Eluc_v4_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.fna.gz

gunzip GCF_004634155.1_Eluc_v4_genomic.fna.gz 
bgzip GCF_004634155.1_Eluc_v4_genomic.fna
samtools faidx GCF_004634155.1_Eluc_v4_genomic.fna.gz

zcat GCF_000233375.1_ICSASG_v2_rm.out.gz | ~/sal_enhancers/homeoblock_alignments/rm_out2bed.py | sort -k1,1 -k2,2n | bgzip -c > salsal_rep_coords.bed.gz
tabix -pbed salsal_rep_coords.bed.gz

gunzip GCF_000233375.1_ICSASG_v2_genomic.fna.gz
bedtools maskfasta -fi GCF_000233375.1_ICSASG_v2_genomic.fna -bed salsal_rep_coords.bed.gz -fo GCF_000233375.1_ICSASG_v2_genomic.sm.fna -soft

bgzip GCF_000233375.1_ICSASG_v2_genomic.sm.fna
samtools faidx GCF_000233375.1_ICSASG_v2_genomic.sm.fna.gz 
```

## Generating input files

```shell script
cut -d ',' -f 3,6 lien_et_al_2016_tableS6.csv | sort -u > chromo_comparisons.csv
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.fna.gz -pre salmon
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_004634155.1_Eluc_v4_genomic.fna.gz -pre pike
```

## Pairwise alignments per block

Pairwise alignments were performed per salmon chromosome, against the whole pike genome.

```shell script
python do_all_pairwise.py -chr_list ref_seq_chromo_key.csv -ref_fa /scratch/project_2002047/sal_enhance/alignments/GCF_004634155.1_Eluc_v4_genomic.rename.fa -ref_name pike -query_fa /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.rename.fa -query_name salmon -out /scratch/project_2002047/sal_enhance/alignments/
ls /scratch/project_2002047/sal_enhance/alignments/NC_0273*/*maf | cut -d '.' -f 1,2,3,4,5 | while read i; do single_cov2 $i.maf R=pike S=pike > $i.sing.maf; done
ls /scratch/project_2002047/sal_enhance/alignments/NC_0273*/*.sing.maf | python make_b_mafs.py 
```

## Multiple alignments per block

Multiple alignments were then generated for each unique pair of salmon chromosomes listed in table S6 in Lien et al (2016) against the pike genome.



