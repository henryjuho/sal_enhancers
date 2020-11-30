# Alignments of salmon homeoblocks with northern pike outgroup (went unused)

## Genomes

Northern pike and atlantic salmon genomes were downloaded from NCBI. Atlantic salmon was then soft masked for repeat regions (pike already masked). Additionally the salmon CDS regions were hard masked to leave only non-coding bases.

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
bedtools maskfasta -fi GCF_000233375.1_ICSASG_v2_genomic.sm.fna -bed GCF_000233375.1_ICSASG_v2_cds.bed.gz -fo GCF_000233375.1_ICSASG_v2_genomic.sm.nc.fna

bgzip GCF_000233375.1_ICSASG_v2_genomic.sm.nc.fna
samtools faidx GCF_000233375.1_ICSASG_v2_genomic.sm.nc.fna.gz 
```

## Generating input files

```shell script
cp /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.nc.fna.gz /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.nc.salb.fna.gz
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.nc.fna.gz -pre salmon
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.nc.salb.fna.gz -pre salmon_b
python fasta_add_header_prefix.py -fa /scratch/project_2002047/sal_enhance/alignments/GCF_004634155.1_Eluc_v4_genomic.fna.gz -pre pike
```

## Pairwise alignments per block

Pairwise alignments were performed per salmon chromosome, against the whole pike genome.

```shell script
python do_all_pairwise.py -block_info lien_et_al_2016_tableS6.csv -key ref_seq_chromo_key.csv -sal_ref /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.nc.rename.fa -sal_query /scratch/project_2002047/sal_enhance/alignments/GCF_000233375.1_ICSASG_v2_genomic.sm.nc.salb.rename.fa -pike_query /scratch/project_2002047/sal_enhance/alignments/GCF_004634155.1_Eluc_v4_genomic.rename.fa -out /scratch/project_2002047/sal_enhance/alignments/
ls /scratch/project_2002047/sal_enhance/alignments/*-*/*maf | cut -d '.' -f 1-2 | while read i; do single_cov2 $i.maf R=salmon S=salmon > $i.sing.maf; done
```

## Multiple alignments per block

Multiple alignments were then generated for each homeoblock listed in table S6 in Lien et al (2016) with the pike genome as outgroup.

```shell script
ls -d /scratch/project_2002047/sal_enhance/alignments/*-* | python chromo_pair_multi.py 
mkdir /scratch/project_2002047/sal_enhance/alignments/block_multiple
cp /scratch/project_2002047/sal_enhance/alignments/*-*/aligned/*maf /scratch/project_2002047/sal_enhance/alignments/block_multiple
```

These were then converted to ```.fasta``` format:

```shell script
python get_block_fastas.py -maf_dir /scratch/project_2002047/sal_enhance/alignments/block_multiple/
python clean_fastas.py -fa_dir /scratch/project_2002047/sal_enhance/alignments/block_multiple/ > homeoblock_alignment_summary.csv
```

## Divergence analysis

Within block divergence was estimated using the K80 model in ape.

```shell script
echo 'block,pairwise,salmon_branch,salmonb_branch' > homeoblock_divergence.csv
ls /scratch/project_2002047/sal_enhance/alignments/block_multiple/*.clean.fa | while read i; do Rscript k80_div_est.R $i; done >> homeoblock_divergence.csv
```

Divergence results: [homeoblock_divergence.csv](homeoblock_divergence.csv)