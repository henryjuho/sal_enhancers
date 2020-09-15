#!/usr/bin/env python

import argparse
from qsub import q_sub
import sys


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', help='Reference genome', required=True)
    parser.add_argument('-out', help='Output dir', required=True)
    args = parser.parse_args()

    for vcf in sys.stdin:

        vcf = vcf.rstrip()
        out_vcf = args.out + vcf.replace('.allsites.g.vcf', '.raw.snps.indels.vcf').split('/')[-1]

        # submit job
        genotyper = ('gatk --java-options "-Xmx4g -Djava.io.tmpdir=/scratch/project_2002047/tmp" GenotypeGVCFs '
                     '-R {} -V {} -O {} ').format(args.ref, vcf, out_vcf)

        q_sub([genotyper], out=out_vcf.replace('.vcf', ''), t=60, mem=10, rmem=10, scheduler='SLURM')


if __name__ == '__main__':
    main()
