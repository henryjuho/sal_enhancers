#!/usr/bin/env python

import argparse
import gzip

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fa', help='Fasta file to add prefix to headers', required=True)
parser.add_argument('-pre', help='Prefix', required=True)
args = parser.parse_args()

# variables
fa = args.fa
prefix = args.pre
out = fa.replace('.fna.gz', '.rename.fa')

# rename
with open(out, 'w') as new_fa:
    for line in gzip.open(fa, 'r'):
        line = line.decode("utf-8")
        if line.startswith('>'):
            line = line.split()[0].replace('>', '>' + prefix + '.') + '\n'
            new_fa.write(line)
        else:
            new_fa.write(line)

print('Done')
