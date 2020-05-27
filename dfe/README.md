# Getting per gene sfs data

For 4fold sites, CDS, UTRs, introns.

```shell script
mkdir /scratch/project_2002047/sal_enhance/sfs

python chromosomal_bootstrap_sfs_data.py -bed_regs /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -bed_target /scratch/tuyida/bartonhe/sal_ref/salmo_salar_4fold.bed.gz -vcf /scratch/project_2002047/sal_reseq/post_vqsr_barson/salsal_30.autosomes.t99_5_snps.allfilters.polarised.phased.vcf.gz -call_fa /scratch/project_2002047/sal_reseq/callable_sites_pol_barson/salsal_30.callable.fa -region 4fold -chromo_list /scratch/tuyida/bartonhe/sal_ref/autosomes_list.txt -out_dir /scratch/project_2002047/sal_enhance/sfs/ 
head -n 1 /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.NC_027300.1.txt > /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt
cat /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.NC*.1.txt | grep -v ^region >> /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt

python chromosomal_bootstrap_sfs_data.py -bed_regs /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -bed_target /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_cds.bed.gz -vcf /scratch/project_2002047/sal_reseq/post_vqsr_barson/salsal_30.autosomes.t99_5_snps.allfilters.polarised.phased.vcf.gz -call_fa /scratch/project_2002047/sal_reseq/callable_sites_pol_barson/salsal_30.callable.fa -region cds -chromo_list /scratch/tuyida/bartonhe/sal_ref/autosomes_list.txt -out_dir /scratch/project_2002047/sal_enhance/sfs/ 
head -n 1 /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.cds.NC_027300.1.txt > /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.cds.all.txt
cat /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.cds.NC_0273*1.txt | grep -v ^region >> /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.cds.all.txt

python chromosomal_bootstrap_sfs_data.py -bed_regs /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -bed_target /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_utrs.bed.gz -vcf /scratch/project_2002047/sal_reseq/post_vqsr_barson/salsal_30.autosomes.t99_5_snps.allfilters.polarised.phased.vcf.gz -call_fa /scratch/project_2002047/sal_reseq/callable_sites_pol_barson/salsal_30.callable.fa -region utr -chromo_list /scratch/tuyida/bartonhe/sal_ref/autosomes_list.txt -out_dir /scratch/project_2002047/sal_enhance/sfs/ 
head -n 1 /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.utr.NC_027300.1.txt > /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.utr.all.txt
cat /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.utr.NC_0273*1.txt | grep -v ^region >> /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.utr.all.txt

python chromosomal_bootstrap_sfs_data.py -bed_regs /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_gene_names.bed.gz -bed_target /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_introns.bed.gz -vcf /scratch/project_2002047/sal_reseq/post_vqsr_barson/salsal_30.autosomes.t99_5_snps.allfilters.polarised.phased.vcf.gz -call_fa /scratch/project_2002047/sal_reseq/callable_sites_pol_barson/salsal_30.callable.fa -region intron -chromo_list /scratch/tuyida/bartonhe/sal_ref/autosomes_list.txt -out_dir /scratch/project_2002047/sal_enhance/sfs/ 
head -n 1 /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.intron.NC_027300.1.txt > /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.intron.all.txt
cat /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.intron.NC_0273*1.txt | grep -v ^region >> /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.intron.all.txt
```

For enhancers, the nearest gene was identified for each region (for bootstrapping by gene).

```shell script
cat filtered_outputPeaks.bed | sort -k1,1 -k2,2n | cut -f1,2,3 | bgzip -c > enhancer_peaks.bed.gz
tabix -pbed enhancer_peaks.bed.gz 

bedtools closest -t first -a enhancer_peaks.bed.gz -b /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_gene_names.bed.gz | cut -f 1,2,3,7 | bgzip -c > enhancer_peaks_nearestgene.bed.gz
tabix -pbed enhancer_peaks_nearestgene.bed.gz
```

Then enhancers were subset into UTRs, intronic and intergenic:

```shell script
bedtools intersect -a enhancer_peaks_nearestgene.bed.gz -b /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_utrs.bed.gz | sort -k1,1 -k2,2n | bgzip -c > enhancers_utr.bed.gz
bedtools intersect -a enhancer_peaks_nearestgene.bed.gz -b /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_introns.bed.gz | sort -k1,1 -k2,2n | bgzip -c > enhancers_intron.bed.gz
bedtools intersect -a enhancer_peaks_nearestgene.bed.gz -b /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_intergenic.bed.gz | sort -k1,1 -k2,2n | bgzip -c > enhancers_intergenic.bed.gz

ls enhancers_*.bed.gz | while read i; do tabix -pbed $i; done
```

SFS data for enhancers was then generated:

```shell script
mkdir /scratch/project_2002047/sal_enhance/enhance_sfs

python chromosomal_bootstrap_sfs_data.py -bed_regs NA -bed_target enhancers_utr.bed.gz -vcf /scratch/project_2002047/sal_reseq/post_vqsr_barson/salsal_30.autosomes.t99_5_snps.allfilters.polarised.phased.vcf.gz -call_fa /scratch/project_2002047/sal_reseq/callable_sites_pol_barson/salsal_30.callable.fa -region utr_enhancers -chromo_list /scratch/tuyida/bartonhe/sal_ref/autosomes_list.txt -out_dir /scratch/project_2002047/sal_enhance/enhance_sfs/
head -n 1 /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.utr_enhancers.NC_027300.1.txt > /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.utr_enhancers.all.txt
cat /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.utr_enhancers.NC*.1.txt | grep -v ^region >> /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.utr_enhancers.all.txt

python chromosomal_bootstrap_sfs_data.py -bed_regs NA -bed_target enhancers_intron.bed.gz -vcf /scratch/project_2002047/sal_reseq/post_vqsr_barson/salsal_30.autosomes.t99_5_snps.allfilters.polarised.phased.vcf.gz -call_fa /scratch/project_2002047/sal_reseq/callable_sites_pol_barson/salsal_30.callable.fa -region utr_intron -chromo_list /scratch/tuyida/bartonhe/sal_ref/autosomes_list.txt -out_dir /scratch/project_2002047/sal_enhance/enhance_sfs/
head -n 1 /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.utr_intron.NC_027300.1.txt > /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intron_enhancers.all.txt
cat /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.utr_intron.NC*.1.txt | grep -v ^region | sed 's/utr_intron/intron_enhancers/' >> /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intron_enhancers.all.txt

python chromosomal_bootstrap_sfs_data.py -bed_regs NA -bed_target enhancers_intergenic.bed.gz -vcf /scratch/project_2002047/sal_reseq/post_vqsr_barson/salsal_30.autosomes.t99_5_snps.allfilters.polarised.phased.vcf.gz -call_fa /scratch/project_2002047/sal_reseq/callable_sites_pol_barson/salsal_30.callable.fa -region intergenic_enhancers -chromo_list /scratch/tuyida/bartonhe/sal_ref/autosomes_list.txt -out_dir /scratch/project_2002047/sal_enhance/enhance_sfs/
head -n 1 /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intergenic_enhancers.NC_027300.1.txt > /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intergenic_enhancers.all.txt
cat /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intergenic_enhancers.NC*.1.txt | grep -v ^region >> /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intergenic_enhancers.all.txt
```

# Estimating the DFE from the unfolded SNP SFS

I fitted a model with a gamma distributed DFE to the intronic, utr and cds uSFS with 4fold SNPs as reference. Mutation 
rates were assumed to be equal between neutral and selected sites. Runs were bootstrapped 100 times, through resampling 
and replacement by gene.

SFS data was prepared:

```shell script
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.cds.all.txt -bs_rep 100 > cds_dfe_data.txt
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.utr.all.txt -bs_rep 100 > utr_dfe_data.txt
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.intron.all.txt -bs_rep 100 > intron_dfe_data.txt
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.utr_enhancers.all.txt -bs_rep 100 > utr_enhancers_dfe_data.txt
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intron_enhancers.all.txt -bs_rep 100 > intron_enhancers_dfe_data.txt
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/enhance_sfs/sfs_ncall_regional.intergenic_enhancers.all.txt -bs_rep 100 > intergenic_enhancers_dfe_data.txt
```

Anavar was run:

```shell script
mkdir /scratch/project_2002047/sal_enhance/cds_dfe
cat cds_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/cds_dfe/ss_cds_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/cds_dfe/*results.txt | python gather_bs_reps.py > salsal30_cds_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/utr_dfe
cat utr_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/utr_dfe/ss_utr_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/utr_dfe/*results.txt | python gather_bs_reps.py > salsal30_utr_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/intron_dfe 
cat intron_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intron_dfe/ss_intron_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intron_dfe/*results.txt | python gather_bs_reps.py > salsal30_intron_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/utr_enhancers_dfe 
cat utr_enhancers_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/ss_utr-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal30_utr-enhancers_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/intron_enhancers_dfe 
cat intron_enhancers_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intron_enhancers_dfe/ss_intron-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intron_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal30_intron-enhancers_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe 
cat intergenic_enhancers_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe/ss_intergenic-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal30_intergenic-enhancers_gamma-dfe_100bs.csv
```

Estimated DFEs were binned into selective categories and 95% confidence intervals calculated:

```shell script
ls salsal30_*_gamma-dfe_100bs.csv | python bin_dfe.py > binned_dfe_allregions.csv
Rscript summarise_dfe.R
```

<img src="all_regions_dfe.png" width="300">
