#!/usr/bin/env python

import argparse
from qsub import *


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', help='Reference genome', required=True)
    parser.add_argument('-bam_dir', help='Directory containing BAM files', action='append', required=True)
    parser.add_argument('-out_dir', help='Output directory', required=True)
    args = parser.parse_args()

    # variables
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    # per bam jobs
    for bam_dir in args.bam_dir:
        bams = [x for x in os.listdir(bam_dir) if x.endswith('.bam')]

        for b in bams:
            sample = b.replace('_dupremRG.bam', '').replace('_merged', '')
            bam = bam_dir + b
            hap_caller = ('gatk --java-options "-Xmx4g" HaplotypeCaller '
                          '-R ' + args.ref + ' '
                          '-I ' + bam + ' '
                          '-ERC GVCF '
                          '-ploidy 2 '
                          '-O ' + args.out_dir + sample + '.raw.snps.indels.g.vcf')

            q_sub([hap_caller], out=args.out_dir + sample + '.hap_calling', t=60, rmem=8, mem=8, scheduler='SLURM')


if __name__ == '__main__':
    main()
