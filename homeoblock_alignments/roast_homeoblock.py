#!/usr/bin/env python

import argparse
import os
import subprocess


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-maf_dir', help='top level directory of pairwise mafs', required=True)
    parser.add_argument('-ref', help='Reference species name', required=True)
    args = parser.parse_args()

    maf_dir = args.maf_dir
    mafs = maf_dir + '*.maf'

    out_dir = maf_dir + 'aligned/'
    out_maf = out_dir + maf_dir.split('/')[-2] + '.multiple.maf'

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    temp_dir = out_dir + 'tmp/'
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    tree = ('"((salmon salmon_b) pike)"')

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
