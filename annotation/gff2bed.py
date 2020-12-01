#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-gene', help='if specified will out put gene name too', default=False, action='store_true')
    args = parser.parse_args()

    for line in sys.stdin:
        if line.startswith('#'):
            continue
        else:
            line = line.split()
            chromo, start, end = line[0], int(line[3])-1, line[4]

            if args.gene:
                info = line[8]
                gene_name = info.split(';')[2].split('=')[1]

                print(chromo, start, end, gene_name, sep='\t')

            else:
                print(chromo, start, end, sep='\t')


if __name__ == '__main__':
    main()
