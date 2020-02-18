#!/usr/bin/env python

import sys


def dict_print(dictionary):

    for k in sorted(dictionary.keys()):
        print(k, ': ', dictionary[k], sep='')

    print()


def main():

    seqs = {'salmon': '', 'salmon_b': '', 'pike': ''}

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
                    missing = ''.join(['M' for i in range(0, seq_len)])
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
        miss_len = seqs[spp].count('M')
        print(spp, align_len, miss_len, round(miss_len/align_len, 3), sep=',')

    # out fasta
    for spp in ('salmon', 'salmon_b', 'pike'):
        print('>' + spp)

        for i in range(0, align_len, 60):
            print(seqs[spp][i: i+60])

        print()


if __name__ == '__main__':
    main()
