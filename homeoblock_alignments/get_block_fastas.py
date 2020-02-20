#!/usr/bin/env python

import argparse
import os
from qsub import q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-maf_dir', help='directory of block maf files', required=True)
    args = parser.parse_args()

    mafs = [x for x in os.listdir(args.maf_dir) if x.endswith('.maf')]

    for maf in mafs:

        block = maf.replace('.multiple.maf', '')
        maf = args.maf_dir + maf

        cmd = ('cat {maf} | python /users/bartonhe/sal_enhancers/homeoblock_alignments/maf2fasta.py '
               '-block {block} -out {out}').format(maf=maf, block=block, out=args.maf_dir)

        q_sub([cmd], out=args.maf_dir + block + '.clean_align', scheduler='SLURM')


if __name__ == '__main__':
    main()
