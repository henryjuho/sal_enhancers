# Estimating the DFE from the unfolded SNP SFS

I fitted a model with a gamma distributed DFE to the intronic, utr and cds uSFS with 4fold SNPs as reference. Mutation 
rates were assumed to be equal between neutral and selected sites. Runs were bootstrapped 100 times, through resampling 
and replacement by gene.


Anavar was run:

```shell script
mkdir /scratch/project_2002047/sal_enhance/cds_dfe
cat ../sfs/cds_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/cds_dfe/ss_cds_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/cds_dfe/*results.txt | python gather_bs_reps.py > salsal31_cds_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/utr_dfe
cat ../sfs/utr_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/utr_dfe/ss_utr_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/utr_dfe/*results.txt | python gather_bs_reps.py > salsal31_utr_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/intron_dfe 
cat ../sfs/intron_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intron_dfe/ss_intron_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intron_dfe/*results.txt | python gather_bs_reps.py > salsal31_intron_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/utr_enhancers_dfe 
cat ../sfs/utr_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/ss_utr-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
# 4 jobs failed: todo
ll /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/*error | grep -vw 0 | tr -s ' ' | cut -d ' ' -f 9 | cut -d '.' -f 1-3 | while read i; do sbatch $i.sh; done
ll /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/*error | grep -w 0 | tr -s ' ' | cut -d ' ' -f 9 | cut -d '.' -f 1-3 | while read i; do echo $i.results.txt; done | python gather_bs_reps.py > salsal31_utr-enhancers_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/intron_enhancers_dfe 
cat ../sfs/intron_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intron_enhancers_dfe/ss_intron-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intron_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal31_intron-enhancers_gamma-dfe_100bs.csv

mkdir /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe 
cat ../sfs/intergenic_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe/ss_intergenic-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal31_intergenic-enhancers_gamma-dfe_100bs.csv
```

Estimated DFEs were binned into selective categories and 95% confidence intervals calculated:

```shell script
ls salsal31_*_gamma-dfe_100bs.csv | python bin_dfe.py > binned_dfe_allregions.csv
Rscript summarise_dfe.R
```

<img src="all_regions_dfe.png" width="300">
