#!/usr/bin/env python

import argparse
from qsub import q_sub
import sys
import os


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', help='Reference genome', required=True)
    parser.add_argument('-out_dir', help='Output directory', required=True)
    args = parser.parse_args()

    chromo_list = [x.split('\t')[0] for x in open(args.ref + '.fai') if not x.startswith('NW')]
    females = ['Uts_11_53', 'Uts_11_52', 'Uts_11_39', 'Uts_11_29', 'Uts_11_28',
               'Naus_12_0037', 'Nams_12_0071', 'Jols_13_0001', 'Arga_12_0082', 'Alta_12_0124']

    # per bam jobs
    bams = [x.rstrip() for x in sys.stdin]

    for b in bams:

        job_list = []
        sample = b.split('/')[-1].split('.')[0]

        # loop through chromos
        for chromo in chromo_list:

            # set ploidy for mito and sdy
            if chromo == 'NC_001960.1' or chromo.startswith('KT'):
                ploidy = 1
            else:
                ploidy = 2

            chromo_dir = args.out_dir + chromo + '/'

            # create chromo dir if not already there
            if not os.path.isdir(chromo_dir):
                os.makedirs(chromo_dir)

            # skip SDY for relevant females
            if sample in females and chromo.startswith('KT'):
                continue

            out_gvcf = chromo_dir + sample + '.' + chromo + '.allsites.g.vcf'

            hap_caller = ('gatk --java-options "-Xmx4g" HaplotypeCaller '
                          '-R {ref} '
                          '-I {bam} '
                          '-ERC GVCF '
                          '-ploidy {ploidy} '
                          '-O {gvcf} '
                          '-L {chr} ').format(ref=args.ref, bam=b, ploidy=ploidy, gvcf=out_gvcf, chr=chromo)

            job_list.append(hap_caller)

        # submit one job per bam
        out_stem = args.out_dir + sample + '.hap_calling'
        q_sub(job_list, out=out_stem, t=60, rmem=8, mem=8, scheduler='SLURM')


if __name__ == '__main__':
    main()
