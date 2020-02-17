#!/usr/bin/env python

import argparse
import os
import subprocess


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed_dir', help='directory of wga beds', required=True)
    parser.add_argument('-block_info', help='csv file of block coordinates and names', required=True)
    parser.add_argument('-key', help='refseq ID chromo name key', required=True)
    parser.add_argument('-out', help='Output directory', required=True)
    args = parser.parse_args()

    keys = {x.split(',')[1].rstrip(): x.split(',')[0] for x in open(args.key)}

    # construct bed comb dictionary
    bed_dict = {}
    for in_bed in os.listdir(args.bed_dir):

        if not in_bed.endswith('.wga.bed.gz'):
            continue

        chr_a, chr_b = in_bed.split('/')[-1].split('.mult')[0].split('_vs_')

        if chr_a not in bed_dict.keys():
            bed_dict[chr_a] = {}

        bed_dict[chr_a][chr_b] = args.bed_dir + in_bed

    # loop through supplementary table from genome paper
    for line in open(args.block_info):

        if line.startswith('Block'):
            continue

        # block info
        no, block, ssa, start, end, ssa_q, start_q, end_q, strand = line.split(',')

        # convert to ref_seq
        chromo = keys[ssa]
        chromo_q = keys[ssa_q]

        # get required wga bed file
        comb_bed = bed_dict[chromo][chromo_q]

        # bed name
        query_bed_name = args.out + block + '.query_coords.bed'

        # gen bed file for query
        print('writing: ' + query_bed_name)

        with open(query_bed_name, 'w') as q_bed:
            print(chromo_q, start_q, end_q, sep='\t', file=q_bed)

        # output name
        out = args.out + block + '.wga.bed.gz'

        # intersect to get bed per block
        cmd = ('tabix {wga_bed} {chromo}:{start}-{stop} | '
               'python2 ~/WGAbed/non_ref_intersect.py -b {query_bed} -q salmon_b -c {query_chromo} | '
               'bgzip -c > {block_bed_out}'
               '').format(wga_bed=comb_bed, chromo=chromo, start=start, stop=end,
                          query_bed=query_bed_name, query_chromo=chromo_q,
                          block_bed_out=out)

        print('running: ' + cmd + '\n')
        subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
