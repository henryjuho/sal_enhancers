#!/usr/bin/env python

import argparse
import os
import subprocess


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-top_dir', help='top level directory of chromo grouped mafs', required=True)
    parser.add_argument('-ref', help='Reference species name', required=True)
    parser.add_argument('-chr_tag', help='ID for ref chromosome', type=int, required=True)
    parser.add_argument('-n_thou', default=0, type=int)
    args = parser.parse_args()

    # get chromo dir
    all_dirs = [args.top_dir + x + '/' for x in os.listdir(args.top_dir) if x.startswith('BrownTrout')]
    chr_index = (args.n_thou * 1000) + (args.chr_tag - 1)
    maf_dir = all_dirs[chr_index]
    mafs = maf_dir + '*.maf'

    out_dir = maf_dir + 'aligned/'
    out_maf = out_dir + maf_dir.split('/')[-2] + '.multiple.maf'

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    temp_dir = out_dir + 'roast_temp'
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    tree = ('"(((((BrownTrout AtlanticSalmon) (ArcticChar (RainbowTrout (CohoSalmon SockeyeSalmon)))) '
            'DanubeSalmon) EuropeanGrayling) NorthernPike)"')

    ref_name = maf_dir.split('/')[-2].split('.')[0]
    
    os.chdir(maf_dir)

    # construct and submit command line
    roast = ('roast + '
             'T=' + temp_dir + ' '
             'E=' + ref_name + ' ' +
             tree + ' ' +
             mafs + ' ' +
             out_maf)

    subprocess.call(roast, shell=True)


if __name__ == '__main__':
    main()
