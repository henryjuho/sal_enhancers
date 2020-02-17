#!/usr/bin/env python

import argparse
import pysam
import subprocess


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed_dir', help='directory of wga beds', required=True)
    parser.add_argument('-block_info', help='csv file of block coordinates and names', required=True)
    parser.add_argument('-key', help='refseq ID chromo name key', required=True)
    parser.add_argument('-ref', help='reference genome', required=True)
    args = parser.parse_args()

    keys = {x.split(',')[1].rstrip(): x.split(',')[0] for x in open(args.key)}

    ref_fa = pysam.FastaFile(args.ref)

    print('block', 'reference_length', 'align_length', 'align_percent', 'length_outgroup', 'percent_outgroup', 'file',
          sep=',')

    # loop through supplementary table from genome paper
    for line in open(args.block_info):

        if line.startswith('Block'):
            continue

        # block info
        no, block, ssa, start, end, ssa_q, start_q, end_q, strand = line.split(',')

        # convert to ref_seq
        chromo = keys[ssa]

        # count bases in bed
        wga_bed = args.bed_dir + block + '.wga.bed.gz'
        count_cmd = ('')
        n_bases = int(subprocess.Popen(count_cmd, shell=True, stdout=subprocess.PIPE
                                       ).communicate()[0].split('\n')[0])

        # count bases in bed with outgroup
        count_og_cmd = ('')
        og_bases = int(subprocess.Popen(count_og_cmd, shell=True, stdout=subprocess.PIPE
                                        ).communicate()[0].split('\n')[0])

        block_seq = ref_fa.fetch(chromo, int(start), int(end))
        ref_bases = len(block_seq) - block_seq.upper().count('N')

        print(block, ref_bases, n_bases, round(n_bases/ref_bases, 3), og_bases, round(og_bases/ref_bases, 3), wga_bed,
              sep=',')


if __name__ == '__main__':
    main()
