# Getting per gene sfs data

For 4fold sites, CDS, UTRs, introns and enhancers.

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

todo: think about how to bootstrap enhancers.


# Estimating the DFE from the unfolded SNP SFS

I fitted a module with a gamma distributed DFE to the intronic, utr and cds uSFS with 4fold SNPs as reference. Mutation 
rates were assumed to be equal between neutral and selected sites. Runs were bootstrapped 100 times, through resampling 
and replacement by gene.

SFS data was prepared:

```shell script
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.cds.all.txt -bs_rep 100 > cds_dfe_data.txt
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.utr.all.txt -bs_rep 100 > utr_dfe_data.txt
python prep_anavar_data.py -sfs_ref /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.4fold.all.txt -sfs_target /scratch/project_2002047/sal_enhance/sfs/sfs_ncall_regional.intron.all.txt -bs_rep 100 > intron_dfe_data.txt
```

Anavar was run:

```shell script
mkdir /scratch/project_2002047/sal_enhance/cds_dfe
cat cds_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/cds_dfe/ss_cds_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000

mkdir /scratch/project_2002047/sal_enhance/utr_dfe
cat utr_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/utr_dfe/ss_utr_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000

mkdir /scratch/project_2002047/sal_enhance/intron_dfe 
cat intron_dfe_data.txt | python enhancer_dfe.py -n 60 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intron_dfe/ss_intron_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
```

