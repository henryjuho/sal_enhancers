import os
from qsub import q_sub
import sys


def main():

    for line in sys.stdin:

        out_dir = line.rstrip() + '/'

        align_dir = out_dir + 'aligned/'
        tmp_dir = align_dir + 'tmp/'

        print('creating: ' + align_dir)
        os.makedirs(align_dir)

        print('creating: ' + tmp_dir)
        os.makedirs(tmp_dir)

        roast = ('~/sal_enhancers/homeoblock_alignments/roast_homeoblock.py '
                 '-maf_dir {} -ref salmon').format(out_dir)

        q_sub([roast], out=out_dir + 'multiple_align', scheduler='SLURM')

        print()


if __name__ == '__main__':
    main()

