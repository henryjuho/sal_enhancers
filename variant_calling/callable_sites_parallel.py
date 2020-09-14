#!/usr/bin/env python

import argparse
from qsub import q_sub
import math
import sys


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--bed_repeats',
                        help='BED file with repeat regions listed',
                        required=True)
    parser.add_argument('-chr_bed',
                        help='Specifies chromosome to extract callable sites for and coords',
                        default='ALL')
    parser.add_argument('-pol',
                        help='If specified will check if site can be polarised, takes a wga bed file',
                        default='None')
    parser.add_argument('-out_pre',
                        help='output path and prefix',
                        required=True)
    args = parser.parse_args()

    # variables
    vcfs = {x.split('/')[-1].split('.')[1] + '.1': x.rstrip() for x in sys.stdin}
    repeat_bed = args.bed_repeats
    chr_bed = [(x.split()[0], x.split()[1], x.split()[2]) for x in open(args.chr_bed)]
    pol = args.pol
    out_pre = args.out_pre
    file_list = open(out_pre + '_falist.txt', 'w')
    jids = []

    # process each chromosome
    for entry in chr_bed:
        chromo, start, stop = entry[0], int(entry[1]), int(entry[2])
        vcf = vcfs[chromo]

        # min chunk size ~10Mb
        no_chunks = int(math.ceil(stop / 10e6))

        chunks = []
        for i in range(0, no_chunks):
            chunk_size = stop / no_chunks

            if i == no_chunks - 1:
                chunks.append((i, i * chunk_size, stop))

            else:
                chunks.append((i, i * chunk_size, ((i+1) * chunk_size)))

        for part in chunks:

            jid = '{}_part{}.sh'.format(chromo, part[0])
            jids.append(jid)

            fa_out = '{}_{}_part{}.fa'.format(out_pre, chromo, part[0])
            print(fa_out, file=file_list)

            cmd_line = ('python ~/sal_enhancers/variant_calling/callable_sites_from_vcf_regional.py '
                        '-vcf {vcf} '
                        '-bed {bed} '
                        '-DF 2 -mean_depth 8 -N 31 '
                        '-chr {ch} '
                        '-start {st_pos} '
                        '-end {end_pos} '
                        '-pol {wga} '
                        '> {out}'
                        '').format(vcf=vcf,
                                   bed=repeat_bed,
                                   ch=chromo,
                                   st_pos=int(part[1]),
                                   end_pos=int(part[2]),
                                   wga=pol,
                                   out=fa_out)
            q_sub([cmd_line], out=fa_out.replace('.fa', ''), t=48, rmem=10, mem=10, scheduler='SLURM')

    file_list.close()


if __name__ == '__main__':
    main()
