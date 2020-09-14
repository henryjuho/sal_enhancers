#!/usr/bin/env python

import argparse
from qsub import q_sub
import os


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_dir', help='Top level input directory', required=True)
    parser.add_argument('-ref', help='Reference genome', required=True)
    parser.add_argument('-out_dir', help='Output directory', required=True)
    args = parser.parse_args()

    contigs = [x for x in os.listdir(args.in_dir) if os.path.isdir(args.in_dir + x)]
    
    for chromo in contigs:

        chromo_path = args.in_dir + chromo + '/'

        vcf_list = [chromo_path + x for x in os.listdir(chromo_path) if x.endswith('.g.vcf') and 'SRR' not in x]

        out = args.out_dir + 'salsal_{}.{}.allsites.g.vcf'.format(len(vcf_list), chromo)

        # submit job
        combine_cmd = ('gatk --java-options "-Xmx50g -Xms50g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CombineGVCFs '
                       '-R {} '
                       '-O {} ').format(args.ref, out)

        for v in vcf_list:
            combine_cmd += '--variant {} '.format(v)

        q_sub([combine_cmd], out=out.replace('.g.vcf', ''), t=60, mem=53, rmem=53, scheduler='SLURM')


if __name__ == '__main__':
    main()
