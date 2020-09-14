#!/usr/bin/env python

import argparse
import subprocess


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-maf', help='MAF file containg alignment, must be compressed', required=True)
    parser.add_argument('-ref_sp', help='Species to use for coordinates in output bed', required=True)
    parser.add_argument('-ref_sizes', help='Sizes file to extract chromosomes for ref species', required=True)
    parser.add_argument('-out', help='Output directory', required=True)
    parser.add_argument('-chr_tag', help='ID for ref chromosome', type=int, required=True)
    parser.add_argument('-n_thou', default=0, type=int)
    args = parser.parse_args()

    # get target chromo
    chromo_list = [x.split()[0] for x in open(args.ref_sizes)]
    chr_index = (args.n_thou * 1000) + (args.chr_tag - 1)
    chromo = chromo_list[chr_index].replace(args.ref_sp + '.', '')

    # construct conversion cmd
    chromo_bed_out = args.out + 'AtlanticSalmon.' + chromo + '.all.wga.bed.gz'
    maf_to_bed_cmd = ('python2 ~/WGAbed/maf_to_bed.py '
                      '-i ' + args.maf + ' '
                      '-r ' + args.ref_sp + ' '
                      '-c ' + chromo + ' | '
                      'sort -T ' + args.out + ' -k1,1 -k2,2n | '
                      'bgzip -c > ' + chromo_bed_out)

    subprocess.call(maf_to_bed_cmd, shell=True)


if __name__ == '__main__':
    main()
