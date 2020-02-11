import argparse
import os


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-chromo_pairs')
    parser.add_argument('-key')
    parser.add_argument('-in_dir')
    args = parser.parse_args()

    keys = {x.split(',')[1].rstrip(): x.split(',')[0] for x in open(args.key)}

    a_maf = 'pike.salmon.sing.maf'
    b_maf = 'pike.salmon_b.sing.maf'

    for line in open(args.chromo_pairs).readlines()[1:]:

        chr_x, chr_y = line.rstrip().split(',')
        chr_x = keys[chr_x]
        chr_y = keys[chr_y]

        maf_1_dir = args.in_dir + '/' + chr_x + '/'
        maf_1 = [x for x in os.listdir(maf_1_dir) if x.endswith('1.sing.maf')][0]

        maf_2_dir = args.in_dir + '/' + chr_y + '/'
        maf_2 = [x for x in os.listdir(maf_2_dir) if x.endswith('b.sing.maf')][0]

        out_dir = args.in_dir + chr_x + '_vs_' + chr_y + '/'

        align_dir = out_dir + 'aligned/'
        tmp_dir = align_dir + 'tmp/'

        out_maf = chr_x + '_vs_' + chr_y + '.maf'

        print('creating: ' + out_dir)
        print('copying: ' + maf_1_dir + maf_1 + ' -> ' + out_dir + a_maf)
        print('copying: ' + maf_2_dir + maf_2 + ' -> ' + out_dir + b_maf)
        print('creating: ' + align_dir)
        print('creating: ' + tmp_dir)

        print('out maf: ' + out_maf)
        print()


if __name__ == '__main__':
    main()

