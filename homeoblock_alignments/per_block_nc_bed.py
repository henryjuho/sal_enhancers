#!/usr/bin/env python

import argparse
import os
import subprocess


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed_dir', help='directory of wga beds', required=True)
    parser.add_argument('-nc_bed', help='bed file of non-coding regions', required=True)
    parser.add_argument('-block_info', help='csv file of block coordinates and names', required=True)
    parser.add_argument('-key', help='refseq ID chromo name key', required=True)
    args = parser.parse_args()

    for line in open(args.block_info):


        cmd = ('')