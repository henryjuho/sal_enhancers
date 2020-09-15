#!/usr/bin/env python

import argparse
from qsub import q_sub
import sys


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', help='Reference genome', required=True)
    parser.add_argument('-out', help='Output vcf stem', required=True)
    args = parser.parse_args()

    # vcfs
    vcf_list = [x.rstrip() for x in sys.stdin]

    contigs = [x.split('\t')[0] for x in open(args.ref + '.fai') if not x.startswith('NW')
               and not x.startswith('KT') and not x.startswith('NC_001960.1')]

    for chromo in contigs:

        out = args.out + '_' + chromo + '.gatk.allsites.g.vcf'

        # submit job
        combine_cmd = ('gatk --java-options "-Xmx20g -Xms20g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CombineGVCFs '
                       '-R {} '
                       '-O {} '
                       '-L {} ').format(args.ref, out, chromo)

        for v in vcf_list:
            combine_cmd += '--variant {} '.format(v)

        q_sub([combine_cmd], out=out.replace('.g.vcf', ''), t=60, mem=25, rmem=25, scheduler='SLURM')


if __name__ == '__main__':
    main()
