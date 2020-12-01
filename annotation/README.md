# Preparation of reference files

## Ref files

RefSeq directory: <ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2>

```bash
mkdir /scratch/tuyida/bartonhe/sal_ref/
cd /scratch/tuyida/bartonhe/sal_ref/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.gff.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_rna.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_cds_from_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_rm.out.gz
```

## Degeneracy annotation

```bash
python annotate_degen.py -cds_fa /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_cds_from_genomic.fna.gz
ls /scratch/tuyida/bartonhe/sal_ref/*bed.gz | while read i; do tabix -pbed $i; done
```

## Genomic regions from gff file

```bash
cd /scratch/tuyida/bartonhe/sal_ref

cat /scratch/project_2002047/sal_reseq/bams/Reference_genome_with_SDY/GCF_000233375.1_ICSASG_v2_genomic_with_SDY.fna.fai | cut -f 1,2 > GCF_000233375.1_ICSASG_v2.genome.txt

zgrep -w exon GCF_000233375.1_ICSASG_v2_genomic.gff.gz | ~/sal_bal_sel/annotation/gff2bed.py | sort -k1,1 -k2,2n | bedtools merge | bgzip -c > GCF_000233375.1_ICSASG_v2_exons.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_exons.bed.gz 

zgrep -Pw 'Gnomon\tgene' GCF_000233375.1_ICSASG_v2_genomic.gff.gz | ~/sal_bal_sel/annotation/gff2bed.py | sort -k1,1 -k2,2n | bedtools merge | bgzip -c > GCF_000233375.1_ICSASG_v2_genes.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_genes.bed.gz 

zgrep -Pw 'Gnomon\tgene' GCF_000233375.1_ICSASG_v2_genomic.gff.gz | ~/sal_bal_sel/annotation/gff2bed.py -gene | sort -k1,1 -k2,2n | bgzip -c > GCF_000233375.1_ICSASG_v2_gene_names.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_gene_names.bed.gz

zgrep -Pw 'Gnomon\tCDS' GCF_000233375.1_ICSASG_v2_genomic.gff.gz | ~/sal_bal_sel/annotation/gff2bed.py | sort -k1,1 -k2,2n | bedtools merge | bgzip -c > GCF_000233375.1_ICSASG_v2_cds.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_cds.bed.gz 

bedtools subtract -a GCF_000233375.1_ICSASG_v2_genes.bed.gz -b GCF_000233375.1_ICSASG_v2_exons.bed.gz | bgzip -c > GCF_000233375.1_ICSASG_v2_introns.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_introns.bed.gz 

bedtools subtract -a GCF_000233375.1_ICSASG_v2_exons.bed.gz -b GCF_000233375.1_ICSASG_v2_cds.bed.gz | bgzip -c > GCF_000233375.1_ICSASG_v2_utrs.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_utrs.bed.gz

bedtools complement -i GCF_000233375.1_ICSASG_v2_genes.bed.gz -g GCF_000233375.1_ICSASG_v2.genome.txt | bgzip -c > GCF_000233375.1_ICSASG_v2_intergenic.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_intergenic.bed.gz 

grep ^NC_02 GCF_000233375.1_ICSASG_v2_genomic.fna.fai |cut -f1,2 | awk 'BEGIN { FS ="\t"}; {print $1, 0, $2}' > GCF_000233375.1_ICSASG_v2_autosomes.bed
```

<!---

These were also produced containing gene name information as follows:

```bash
bedtools intersect -a GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -b salmo_salar_0fold.bed.gz | sort -k1,1 -k2,2n | bgzip -c > salmo_salar_0fold_genenames.bed.gz 
tabix -pbed salmo_salar_0fold_genenames.bed.gz 
bedtools intersect -a GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -b salmo_salar_4fold.bed.gz | sort -k1,1 -k2,2n | bgzip -c > salmo_salar_4fold_genenames.bed.gz 
tabix -pbed salmo_salar_4fold_genenames.bed.gz 
bedtools intersect -a GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -b GCF_000233375.1_ICSASG_v2_cds.bed.gz | sort -k1,1 -k2,2n | bgzip -c > GCF_000233375.1_ICSASG_v2_cds_genenames.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_cds_genenames.bed.gz 
bedtools intersect -a GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -b GCF_000233375.1_ICSASG_v2_introns.bed.gz | sort -k1,1 -k2,2n | bgzip -c > GCF_000233375.1_ICSASG_v2_introns_genenames.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_introns_genenames.bed.gz 
```
--->

## Getting a bed file of repeats

```bash
cd /scratch/tuyida/bartonhe/sal_ref/
zcat GCF_000233375.1_ICSASG_v2_rm.out.gz | python /users/bartonhe/sal_bal_sel/annotation/rm_out2bed.py | bgzip -c > salmo_salar_repeats.bed.gz
tabix -p bed salmo_salar_repeats.bed.gz
```

## Duplicated regions

Coordinates of homologous genomic blocks were obtained from supplementary table 6 from Lien et al. (2016). And used to generate a bed file of duplicated regions.

```bash
zgrep -Pw 'RefSeq\tregion' /scratch/tuyida/bartonhe/sal_ref/*gff.gz | tr -s " " | cut -f 1,9 | python -c "import sys; [print(x.split()[0], 'ssa' + x.split('ssa')[1].split(';')[0], sep=',') for x in sys.stdin if x.startswith('NC_02')]" > /scratch/tuyida/bartonhe/sal_ref/ref_seq_chromo_key.csv
cat duplicated_regions.csv | python dup_csv2bed.py  -key /scratch/tuyida/bartonhe/sal_ref/ref_seq_chromo_key.csv | sort -k1,1 -k2,2n | bgzip -c > /scratch/tuyida/bartonhe/sal_ref/duplicated_region_blocks.bed.gz
tabix -pbed /scratch/tuyida/bartonhe/sal_ref/duplicated_region_blocks.bed.gz 
```

These were used to get coordinates for annotated sites not in these regions:

```bash
cd /scratch/tuyida/bartonhe/sal_ref/

bedtools subtract -a GCF_000233375.1_ICSASG_v2_intergenic.bed.gz -b duplicated_region_blocks.bed.gz | bgzip -c > GCF_000233375.1_ICSASG_v2_intergenic_notdup.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_intergenic_notdup.bed.gz 

bedtools subtract -a GCF_000233375.1_ICSASG_v2_introns.bed.gz -b duplicated_region_blocks.bed.gz | bgzip -c > GCF_000233375.1_ICSASG_v2_introns_notdup.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_introns_notdup.bed.gz 

bedtools subtract -a GCF_000233375.1_ICSASG_v2_cds.bed.gz -b duplicated_region_blocks.bed.gz | bgzip -c > GCF_000233375.1_ICSASG_v2_cds_notdup.bed.gz
tabix -pbed GCF_000233375.1_ICSASG_v2_cds_notdup.bed.gz 

bedtools subtract -a salmo_salar_4fold.bed.gz -b duplicated_region_blocks.bed.gz | bgzip -c > salmo_salar_4fold_notdup.bed.gz
tabix -pbed salmo_salar_4fold_notdup.bed.gz

bedtools subtract -a salmo_salar_0fold.bed.gz -b duplicated_region_blocks.bed.gz | bgzip -c > salmo_salar_0fold_notdup.bed.gz
tabix -pbed salmo_salar_0fold_notdup.bed.gz
```
