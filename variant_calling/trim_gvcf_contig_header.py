#!/usr/bin/env python

import argparse
from qsub import q_sub
import os


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_dir', help='Top level input directory', required=True)
    parser.add_argument('-out_dir', help='Output directory', required=True)
    args = parser.parse_args()

    contigs = [x for x in os.listdir(args.in_dir) if os.path.isdir(args.in_dir + x)]

    for chromo in contigs:

        chromo_path = args.in_dir + chromo + '/'
        out_chromo_path = args.out_dir + chromo + '/'

        if not os.path.isdir(out_chromo_path):
            os.makedirs(out_chromo_path)

        vcf_list = [x for x in os.listdir(chromo_path) if x.endswith('.g.vcf')]

        cmds = []
        for vcf in vcf_list:

            trim_cmd = 'grep -v NW_ {} > {}'.format(chromo_path + vcf, out_chromo_path + vcf)
            cmds.append(trim_cmd)

            index = 'gatk IndexFeatureFile -F ' + out_chromo_path + vcf
            cmds.append(index)

        q_sub(cmds, out=args.out_dir + chromo + '_trimhead', t=10, scheduler='SLURM')


if __name__ == '__main__':
    main()
