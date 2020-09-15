# Salmonid whole genome alignment

## Setup and genomes

Genomes for whole genome alignment were downloaded as follows:

```bash
mkdir /scratch/tuyida/bartonhe/sal_alignment
cd /scratch/tuyida/bartonhe/sal_alignment

mkdir genomes
cd genomes/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/233/375/GCF_000233375.1_ICSASG_v2/GCF_000233375.1_ICSASG_v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/163/495/GCF_002163495.1_Omyk_1.0/GCF_002163495.1_Omyk_1.0_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/021/735/GCF_002021735.1_Okis_V1/GCF_002021735.1_Okis_V1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/910/315/GCF_002910315.2_ASM291031v2/GCF_002910315.2_ASM291031v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/901/001/165/GCF_901001165.1_fSalTru1.1/GCF_901001165.1_fSalTru1.1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/317/085/GCA_003317085.1_ASM331708v1/GCA_003317085.1_ASM331708v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/348/285/GCA_004348285.1_ASM434828v1/GCA_004348285.1_ASM434828v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/634/155/GCF_004634155.1_Eluc_v4/GCF_004634155.1_Eluc_v4_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/149/115/GCF_006149115.1_Oner_1.0/GCF_006149115.1_Oner_1.0_genomic.fna.gz

ls *.gz > spp_names.txt
```

Fasta contig headers renamed to include species info at start (```spp_names.txt``` had second column manually added containing common spp names):

## Pairwise alignments

Pairwise alignments where performed using lastz using brown trout as the reference genome, 
due to its better assembly and close position to Atlantic salmon phylogenetically. These alignments 
were performed per ref chromosome. So the whole query genome aligned against each brown trout contig. 

The alignment pipeline follows the UCSC pipeline outlined here: <http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto>.

```bash
python rename_fastas.py -fa_list /scratch/tuyida/bartonhe/sal_alignment/genomes/spp_names.txt
python do_all_pairwise.py -fa_list /scratch/tuyida/bartonhe/sal_alignment/genomes/spp_names.txt -ref_spp BrownTrout -out /scratch/tuyida/bartonhe/sal_alignment/pairwise/
```

## Multiple alignment

```bash
ls -d /scratch/tuyida/bartonhe/sal_alignment/pairwise/BrownTrout*/ | python all_single_cov.py 
python group_pairwise_by_chr.py -top_dir /scratch/tuyida/bartonhe/sal_alignment/pairwise/ -chr_list /scratch/tuyida/bartonhe/sal_alignment/pairwise/BrownTrout.AtlanticSalmon/BrownTrout.sizes -out_dir /scratch/tuyida/bartonhe/sal_alignment/multiple/
python roast_all.py -top_dir /scratch/tuyida/bartonhe/sal_alignment/multiple/ -ref BrownTrout -out /scratch/tuyida/bartonhe/sal_alignment/multiple/
mkdir /scratch/tuyida/bartonhe/sal_alignment/wgabed
ls /scratch/tuyida/bartonhe/sal_alignment/multiple/BrownTrout.N*/aligned/*.maf | python merge_mafs.py -out_maf /scratch/tuyida/bartonhe/sal_alignment/wgabed/BrownTrout.8_query_spp.maf.gz
```

## Converting from MAF to WGAbed

The 29 assembled Atlantic salmon chromosomes (excluding mitochondrion) were extracted from the ```.maf``` file and converted to a ```.wga.bed``` file. 

```bash
mkdir /scratch/tuyida/bartonhe/sal_alignment/wgabed/chromosomal_atlantic_salmon
python into_bed.py -maf /scratch/tuyida/bartonhe/sal_alignment/wgabed/BrownTrout.8_query_spp.maf.gz -ref_sp AtlanticSalmon -ref_sizes /scratch/tuyida/bartonhe/sal_alignment/pairwise/BrownTrout.AtlanticSalmon/AtlanticSalmon.sizes -out /scratch/tuyida/bartonhe/sal_alignment/wgabed/chromosomal_atlantic_salmon/ -no_sub
sbatch --account=tuyida /scratch/tuyida/bartonhe/sal_alignment/wgabed/chromosomal_atlantic_salmon/all_wgabed_job.sh 

cd /scratch/tuyida/bartonhe/sal_alignment/wgabed/chromosomal_atlantic_salmon
zcat AtlanticSalmon.NC_0273*.bed.gz | bgzip -c > AtlanticSalmon.9way.wga.bed.gz
tabix -pbed AtlanticSalmon.9way.wga.bed.gz 
mv AtlanticSalmon.9way.wga.bed.gz* ../
cd ~/sal_bal_sel/genome_alignment/
zcat /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz | python summarise_alignment.py > align_sum.csv
```

Coverage summary [here](align_sum.csv).

## Polarisation

I polarised SNPs using the Atlantic salmon, Arctic charr and brown trout sequences in the alignment and a maximum parsimony approach.

```shell script
python polarise_vcf.py -vcf /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.t99_5_snps.allfilters.vcf.gz -wga_bed /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz

Total no INDELs  : 3723849
INDELs polarised : 1614400
Hotspots         : 632084
Low spp coverage : 1272810
Ambiguous        : 204555
Total unpolarised: 2109449

bgzip /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.t99_5_snps.allfilters.polarised.vcf
tabix -pvcf /scratch/project_2002047/barson_reseq/post_vqsr/salsal_31.autosomes.t99_5_snps.allfilters.polarised.vcf.gz
```