import argparse
from qsub import q_print as q_sub


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-wga', help='whole genome alignment bed', required=True)
    parser.add_argument('-bed', help='bed file of region', required=True)
    parser.add_argument('-region', help='region label', required=True)
    parser.add_argument('-spp', help='Comma separated species list for output', required=True)
    parser.add_argument('-out_stem', help='output directory and file stem', required=True)
    args = parser.parse_args()

    out_stem = args.out_stem + '_' + args.region

    bed_cmd = 'bedtools intersect -a {wga} -b {bed} | bgzip -c > {out_stem}.wga.bed.gz'.format(
        wga=args.wga, bed=args.bed, out_stem=out_stem)

    fasta_cmd = 'zcat {out_stem}.wga.bed.gz | ~/sal_enhancers/divergence/wga2fa.py -out_stem {out_stem}.wga'.format(
        out_stem=out_stem)

    ape_cmd = 'Rscript ~/sal_enhancers/divergence/k80_div_est.R {out_stem}.wga.fa {region} -spp {spp}'.format(
        out_stem=out_stem, region=args.region, spp=args.spp)

    q_sub([bed_cmd, fasta_cmd, ape_cmd], out=args.out_stem, rmem=12, mem=12, scheduler='SLURM')


if __name__ == '__main__':
    main()
