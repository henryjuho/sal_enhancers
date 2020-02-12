import argparse
import os
import shutil
from qsub import q_sub


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-chromo_pairs', help='CSV of chromo comparissions', required=True)
    parser.add_argument('-key', help='CSV of chromosome names and corresponding refseq IDs', required=True)
    parser.add_argument('-in_dir', help='top level input directory', required=True)
    args = parser.parse_args()

    keys = {x.split(',')[1].rstrip(): x.split(',')[0] for x in open(args.key)}

    a_maf = 'pike.salmon.sing.maf'
    b_maf = 'pike.salmon_b.sing.maf'

    for line in open(args.chromo_pairs).readlines()[1:]:

        chr_x, chr_y = line.rstrip().split(',')
        chr_x = keys[chr_x]
        chr_y = keys[chr_y]

        maf_1_dir = args.in_dir + chr_x + '/'
        maf_1 = [x for x in os.listdir(maf_1_dir) if x.endswith('1.sing.maf')][0]

        maf_2_dir = args.in_dir + chr_y + '/'
        maf_2 = [x for x in os.listdir(maf_2_dir) if x.endswith('b.sing.maf')][0]

        out_dir = args.in_dir + chr_x + '_vs_' + chr_y + '/'

        align_dir = out_dir + 'aligned/'
        tmp_dir = align_dir + 'tmp/'

        out_maf = chr_x + '_vs_' + chr_y + '.maf'

        print('creating: ' + out_dir)
        os.makedirs(out_dir)

        print('copying: ' + maf_1_dir + maf_1 + ' -> ' + out_dir + a_maf)
        shutil.copy(maf_1_dir + maf_1, out_dir + a_maf)

        print('copying: ' + maf_2_dir + maf_2 + ' -> ' + out_dir + b_maf)
        shutil.copy(maf_2_dir + maf_2, out_dir + b_maf)

        print('creating: ' + align_dir)
        os.makedirs(align_dir)

        print('creating: ' + tmp_dir)
        os.makedirs(tmp_dir)

        print('out maf: ' + out_dir + out_maf)
        print()

        roast = ('~/sal_enhancers/homeoblock_alignments/roast_homeoblock.py '
                 '-maf_dir {} -ref pike').format(out_dir)

        q_sub([roast], out=out_dir + chr_x + '_vs_' + chr_y, scheduler='SLURM')


if __name__ == '__main__':
    main()

