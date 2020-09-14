#!/usr/bin/env python

import argparse
import sys
import gzip

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-out_maf', help='Output directory and maf name', required=True)
args = parser.parse_args()

maf_list = [x.rstrip() for x in sys.stdin]

# write output
first_maf = True
with gzip.open(args.out_maf, 'wb') as merged_maf:
    for maf in maf_list:
        with open(maf) as maf_in:
            for line in maf_in:
                if line.startswith('#'):
                    if first_maf is True:
                        if line.startswith('##'):
                            merged_maf.write(line.encode('utf-8'))
                else:
                    merged_maf.write(line.encode('utf-8'))
            first_maf = False
            print(maf + ' processed!')

print('Processing complete, merged maf written to ' + args.out_maf)
