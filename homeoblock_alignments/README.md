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
cp /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.fna.gz /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.salb.fna.gz
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.fna.gz -pre salmon
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.salb.fna.gz -pre salmon_b
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_004634155.1_Eluc_v4_genomic.fna.gz -pre pike
```

## Pairwise alignments per block

Pairwise alignments were performed per salmon chromosome, against the whole pike genome.

```shell script
python do_all_pairwise.py -chromo_pairs chromo_comparisons.csv -key ref_seq_chromo_key.csv -sal_ref /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.rename.fa -sal_query /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.salb.rename.fa -pike_query /scratch/project_2002047/sal_enhance/alignments/GCF_004634155.1_Eluc_v4_genomic.rename.fa -out /scratch/project_2002047/sal_enhance/alignments/
ls /scratch/project_2002047/sal_enhance/alignments/*_vs_*/*maf | cut -d '.' -f 1-4 | while read i; do single_cov2 $i.maf R=salmon S=salmon > $i.sing.maf; done
```

## Multiple alignments per block

Multiple alignments were then generated for each unique pair of salmon chromosomes listed in table S6 in Lien et al (2016) against the pike genome.

```shell script
ls -d /scratch/project_2002047/sal_enhance/alignments/NC_0273* | python chromo_pair_multi.py 
updated to here
python chromo_pair_multi.py -chromo_pairs chromo_comparisons.csv -key ref_seq_chromo_key.csv -in_dir /scratch/project_2002047/sal_enhance/alignments/ 
mkdir /scratch/project_2002047/sal_enhance/alignments/all_combinations
cp /scratch/project_2002047/sal_enhance/alignments/*_vs_*/aligned/*maf /scratch/project_2002047/sal_enhance/alignments/all_combinations/
ls /scratch/project_2002047/sal_enhance/alignments/all_combinations/*maf | while read i; do gzip $i; done
```

These were then converted to ```.wga.bed``` format:

```shell script
ls /scratch/project_2002047/sal_enhance/alignments/all_combinations/*maf.gz | python in2bed.py 
ls /scratch/project_2002047/sal_enhance/alignments/all_combinations/*bed.gz | while read i; do tabix -pbed $i; done
```


