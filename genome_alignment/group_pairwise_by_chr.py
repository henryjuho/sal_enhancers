#!/usr/bin/env python

import argparse
import os
import shutil


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-top_dir', help='Directory containing maf files to ensure single coverage for', required=True)
    parser.add_argument('-chr_list', help='ref .sizes file', required=True)
    parser.add_argument('-out_dir', help='Output directory', required=True)
    args = parser.parse_args()

    single_cov_dirs = [args.top_dir + x.rstrip() + '/single_coverage/'
                       for x in os.listdir(args.top_dir) if x.startswith('BrownTrout.')]

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    for line in open(args.chr_list):

        chromo = line.split('\t')[0]
        chromo_dir = args.out_dir + chromo + '/'

        if not os.path.isdir(chromo_dir):
            os.makedirs(chromo_dir)

        for spp_comb in single_cov_dirs:

            spp1, spp2 = spp_comb.split('/')[-3].split('.')

            maf_name = spp_comb + '.'.join([chromo, spp2, 'sing.maf'])
            new_maf = chromo_dir + spp1 + '.' + spp2 + '.sing.maf'

            shutil.copy(maf_name, new_maf)





if __name__ == '__main__':
    main()
