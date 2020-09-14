#!/usr/bin/env python

import argparse
import os
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='Directory containing maf files to ensure single coverage for', required=True)
parser.add_argument('-ref_name', help='Name of reference species', required=True)
args = parser.parse_args()

# variables
directory = args.dir
maf_list = [maf for maf in os.listdir(directory) if maf.endswith('.maf')]
out_dir = directory + 'single_coverage/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
ref = args.ref_name

# run single_cov2
n_mafs = len(maf_list)
counter = 0
for maf in maf_list:
    counter += 1
    in_maf = directory + maf
    output = out_dir + maf.rstrip('.maf') + '.sing.maf'
    cmd_line = 'single_cov2 {maf} [R={ref_spp}] [S={ref_spp}] > {out}'.format(maf=in_maf, ref_spp=ref, out=output)
    subprocess.call(cmd_line, shell=True)
    print(str(counter), '/', str(n_mafs), ' complete')
