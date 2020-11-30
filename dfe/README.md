# Estimating the DFE from the unfolded SNP SFS

I fitted a model with a gamma distributed DFE to the intergenic, intronic, UTR, CDS and 0fold uSFS with 4fold SNPs as 
reference, as well peaks in these regions, and all peaks together. Mutation rates were assumed to be equal between 
neutral and selected sites - a requirement to calculate alpha. Runs were bootstrapped 100 times, through resampling with 
replacement by gene.


## 0 fold degenerate
```shell script
mkdir /scratch/project_2002047/sal_enhance/0fold_dfe
cat ../sfs/0fold_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/0fold_dfe/ss_0fold_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/0fold_dfe/*results.txt | python gather_bs_reps.py > salsal31_0fold_gamma-dfe_100bs.csv
```

## CDS regions

```shell script
mkdir /scratch/project_2002047/sal_enhance/cds_dfe
cat ../sfs/cds_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/cds_dfe/ss_cds_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/cds_dfe/*results.txt | python gather_bs_reps.py > salsal31_cds_gamma-dfe_100bs.csv
```

## UTR

```shell script
mkdir /scratch/project_2002047/sal_enhance/utr_dfe
cat ../sfs/utr_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/utr_dfe/ss_utr_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/utr_dfe/*results.txt | python gather_bs_reps.py > salsal31_utr_gamma-dfe_100bs.csv
```

## Introns

```shell script
mkdir /scratch/project_2002047/sal_enhance/intron_dfe 
cat ../sfs/intron_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intron_dfe/ss_intron_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intron_dfe/*results.txt | python gather_bs_reps.py > salsal31_intron_gamma-dfe_100bs.csv
```

## Intergenic

```shell script
mkdir /scratch/project_2002047/sal_enhance/intergenic_dfe 
cat ../sfs/intergenic_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intergenic_dfe/ss_intergenic_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intergenic_dfe/*results.txt | python gather_bs_reps.py > salsal31_intergenic_gamma-dfe_100bs.csv
```

## Peaks (all)

```shell script
mkdir /scratch/project_2002047/sal_enhance/all_enhancers_dfe
cat ../sfs/all_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/all_enhancers_dfe/ss_all-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/all_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal31_all-enhancers_gamma-dfe_100bs.csv
```

## Peaks (0fold)

```shell script
mkdir /scratch/project_2002047/sal_enhance/0fold_enhancers_dfe  
cat ../sfs/0fold_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/0fold_enhancers_dfe/ss_cds-enhancers_0fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
cd /scratch/project_2002047/sal_enhance/0fold_enhancers_dfe
ls *results.txt | cut -d '_' -f 4-8 | while read i; do mv ss_cds-enhancers_0fold_$i ss_0fold-enhancers_$i; done
cd -
ll /scratch/project_2002047/sal_enhance/0fold_enhancers_dfe/*results.txt | grep -vw 1018 | tr -s ' ' | cut -d ' ' -f 8 | python gather_bs_reps.py > salsal31_0fold-enhancers_gamma-dfe_100bs.csv
```

## Peaks (CDS)

```shell script
mkdir /scratch/project_2002047/sal_enhance/cds_enhancers_dfe  
cat ../sfs/cds_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/cds_enhancers_dfe/ss_cds-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ll /scratch/project_2002047/sal_enhance/cds_enhancers_dfe/*error |  grep -w 0 | tr -s ' ' | cut -d ' ' -f 9 | cut -d '.' -f 1-3 | while read i; do echo $i.results.txt; done | python gather_bs_reps.py > salsal31_cds-enhancers_gamma-dfe_100bs.csv
```

## Peaks (UTR)

```shell script
mkdir /scratch/project_2002047/sal_enhance/utr_enhancers_dfe 
cat ../sfs/utr_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/ss_utr-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ll /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/*error | grep -vw 0 | tr -s ' ' | cut -d ' ' -f 8 | cut -d '.' -f 1-3 | while read i; do sbatch $i.sh; done
ll /scratch/project_2002047/sal_enhance/utr_enhancers_dfe/*error | grep -w 0 | tr -s ' ' | cut -d ' ' -f 8 | cut -d '.' -f 1-3 | while read i; do echo $i.results.txt; done | python gather_bs_reps.py > salsal31_utr-enhancers_gamma-dfe_100bs.csv
```

## Peaks (introns)

```shell script
mkdir /scratch/project_2002047/sal_enhance/intron_enhancers_dfe 
cat ../sfs/intron_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intron_enhancers_dfe/ss_intron-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intron_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal31_intron-enhancers_gamma-dfe_100bs.csv
```

## Peaks (intergenic)

```shell script
mkdir /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe 
cat ../sfs/intergenic_enhancers_sfs_data.txt | python enhancer_dfe.py -n 62 -c 1 -dfe continuous -out_pre /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe/ss_intergenic-enhancers_4fold_continuous_equal_t -constraint equal_mutation_rate -n_search 1000
ls /scratch/project_2002047/sal_enhance/intergenic_enhancers_dfe/*results.txt | python gather_bs_reps.py > salsal31_intergenic-enhancers_gamma-dfe_100bs.csv
```

Estimated DFEs were binned into selective categories and 95% confidence intervals calculated:

```shell script
ls salsal31_*_gamma-dfe_100bs.csv | python bin_dfe.py > binned_gammadfe_allregions_nes.csv
Rscript summarise_dfe.R
```

<img src=binned_gammadfe_allregions_nes.png height=320>