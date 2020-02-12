#!/usr/bin/env python

import sys
from qsub import q_sub


def main():

    # get target chromo
    for line in sys.stdin:

        maf = line.rstrip()
        chromo = maf.split('/')[-1].split('_vs_')[0]
        out_dir = '/'.join(maf.split('/')[:-1])

        # construct conversion cmd
        chromo_bed_out = maf.replace('.maf.gz', '.wga.bed.gz')
        outs = maf.replace('.maf.gz', '.2bed')

        maf_to_bed_cmd = ('python2 ~/WGAbed/maf_to_bed.py '
                          '-i ' + maf + ' '
                          '-r salmon '
                          '-c ' + chromo + ' | '
                          'sort -T ' + out_dir + ' -k1,1 -k2,2n | '
                          'bgzip -c > ' + chromo_bed_out)

        q_sub([maf_to_bed_cmd], out=outs, scheduler='SLURM')


if __name__ == '__main__':
    main()
