#!/usr/bin/env python

import argparse
import sys


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-block', help='homeoblock id', required=True)
    parser.add_argument('-out', help='output directory', required=True)
    args = parser.parse_args()

    seqs = {'salmon': '', 'salmon_b': '', 'pike': ''}

    # out_files
    log = open(args.out + args.block + '.align_summary.csv', 'w')
    fa_out = open(args.out + args.block + '.clean.fa', 'w')

    # process piped maf
    align_block = {}
    first = True

    for line in sys.stdin:
        # skip header
        if line.startswith('#'):
            continue

        # write old block and store new
        elif line.startswith('a'):
            if first:
                first = False
                continue

            # process previous block
            seq_len = max([len(x) for x in align_block.values()])

            for spp in seqs.keys():

                if spp in align_block.keys():
                    seqs[spp] += align_block[spp]

                else:
                    missing = ''.join(['N' for i in range(0, seq_len)])
                    seqs[spp] += missing

            align_block = {}

        elif line.startswith('s'):
            species = line.split()[1].split('.')[0]
            seq = line.rstrip().split()[6]
            align_block[species] = seq

    # process last block
    seq_len = max([len(x) for x in align_block.values()])

    for spp in seqs.keys():

        if spp in align_block.keys():
            seqs[spp] += align_block[spp]

        else:
            missing = ''.join(['N' for i in range(0, seq_len)])
            seqs[spp] += missing

    # summarise
    for spp in ('salmon', 'salmon_b', 'pike'):
        align_len = len(seqs[spp])
        miss_len = seqs[spp].count('N')
        print(args.block, spp, align_len, miss_len, round(miss_len/align_len, 3), sep=',', file=log)

    # clean
    clean_seqs = {'salmon': '', 'salmon_b': '', 'pike': ''}
    for i in range(0, align_len):
        pos_seqs = [seqs['salmon'][i], seqs['salmon_b'][i], seqs['pike']]

        if 'N' in pos_seqs or '-' in pos_seqs:
            continue

        clean_seqs['salmon'] += pos_seqs[0]
        clean_seqs['salmon_b'] += pos_seqs[1]
        clean_seqs['pike'] += pos_seqs[2]

    clean_len = len(clean_seqs['salmon'])
    print(args.block, 'all_cleaned', align_len, clean_len, round(clean_len / align_len, 3), sep=',', file=log)

    # out fasta
    for spp in ('salmon', 'salmon_b', 'pike'):
        print('>' + spp, file=fa_out)

        for i in range(0, align_len, 60):
            print(clean_seqs[spp][i: i+60], file=fa_out)

        print(file=fa_out)

    fa_out.close()
    log.close()


if __name__ == '__main__':
    main()
