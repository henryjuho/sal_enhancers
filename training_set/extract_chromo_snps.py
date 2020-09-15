import sys
from qsub import q_sub
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', help='reference genome', required=True)
    args = parser.parse_args()

    chromo_vcfs = [x.rstrip() for x in sys.stdin]

    for chr_vcf in chromo_vcfs:

        new_vcf = chr_vcf.replace('.allsites', '.raw.snps')

        snp_cmd = ('gatk SelectVariants '
                   '-R {ref} -V {vcf} -O {out} '
                   '--select-type-to-include SNP --exclude-non-variants'
                   '').format(ref=args.ref, vcf=chr_vcf, out=new_vcf)

        q_sub([snp_cmd], out=new_vcf.replace('.vcf', ''), scheduler='SLURM')


if __name__ == '__main__':
    main()


