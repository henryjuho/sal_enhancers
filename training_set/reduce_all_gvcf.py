#!/usr/bin/env python

import argparse
from qsub import q_sub
import sys


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-out', help='Output dir', required=True)
    args = parser.parse_args()

    # vcfs
    vcf_list = [x.rstrip() for x in sys.stdin]

    for vcf in vcf_list:

        new_vcf = args.out + vcf.replace('.g.vcf', '.autosomes.g.vcf').split('/')[-1]

        cmd = 'cat {} | python ~/sal_enhancers/training_set/extract_autosomes.py > {}'.format(vcf, new_vcf)

        q_sub([cmd], out=new_vcf.replace('.g.vcf', ''), rmem=4, mem=4, scheduler='SLURM')


if __name__ == '__main__':
    main()
