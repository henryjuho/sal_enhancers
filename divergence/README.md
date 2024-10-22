# Divergence estimates and alpha

Pipeline for intersecting the whole genome alignment with desired sub regions, converting to a fasta file and 
calculating divergence for the Atlantic salmon branch from common ancestor with brown trout is as follows:

```shell script
mkdir /scratch/project_2002047/sal_enhance/divergence

python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed /scratch/tuyida/bartonhe/sal_ref/salmo_salar_4fold.bed.gz -region 4fold -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_intergenic.bed.gz -region intergenic -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_introns.bed.gz -region intron -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_utrs.bed.gz -region utr -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed /scratch/tuyida/bartonhe/sal_ref/GCF_000233375.1_ICSASG_v2_cds.bed.gz -region cds -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed /scratch/tuyida/bartonhe/sal_ref/salmo_salar_0fold.bed.gz -region 0fold -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div

python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed ~/sal_enhancers/sfs/enhancer_peaks.bed.gz -region all_enhancers -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed ~/sal_enhancers/sfs/enhancers_cds.bed.gz -region cds_enhancers -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed ~/sal_enhancers/sfs/enhancers_utr.bed.gz -region utr_enhancers -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed ~/sal_enhancers/sfs/enhancers_intron.bed.gz -region intron_enhancers -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed ~/sal_enhancers/sfs/enhancers_intergenic.bed.gz -region intergenic_enhancers -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div
python get_region_fasta.py -wga /scratch/tuyida/bartonhe/sal_alignment/wgabed/AtlanticSalmon.9way.wga.bed.gz -bed ~/sal_enhancers/sfs/enhancers_0fold.bed.gz -region 0fold_enhancers -spp AtlanticSalmon,ArcticChar,BrownTrout -out_stem /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div

head -n 1 /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div_all_enhancers.div.txt  > all_regs_div.txt
cat /scratch/project_2002047/sal_enhance/divergence/alpha_calc_div*txt | grep -v region >> all_regs_div.txt
```

Alpha was then calculated using the DFE and the divergence estimates:

```shell script
cat ../dfe/salsal31_* | python calc_alpha_sal.py -div all_regs_div.txt > alpha_estimates.csv
Rscript plot_alpha.R
```

