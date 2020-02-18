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

            dict_print(seqs)

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
            missing = ''.join(['M' for i in range(0, seq_len)])
            seqs[spp] += missing

    # summarise


    # trim out missing (N and M) and indels and output
    #zip()


if __name__ == '__main__':
    main()