#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bam_list', help='List of bam files to genotype', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-out', help='Output_directory and prefix', required=True)
args = parser.parse_args()

# variables
bams = args.bam_list
ref_genome = args.ref
output = args.out + '.samtools.allsites.vcf'

# SAMtools per chromosome

# generate chromosome list
chromosome_list = [x.split('\t')[0] for x in open(args.ref + '.fai')
                   if not x.startswith('NW') and not x.startswith('KT') and not x.startswith('NC_001960.1')]

# generate samtools job for each chromosome
for position in chromosome_list:
    new_output = args.out + '.samtools.allsites.' + position + '.vcf'
    job_name = 'samtools.' + position + '.sh'
    # SAMtools
    SAM_commandline = ('samtools mpileup '
                       '-b ' + bams + ' '
                       '-C 50 '
                       '-f ' + ref_genome + ' '
                       '-r ' + position + ' '
                       '-u '
                       '| bcftools call '
                       '-O v '
                       '-m '
                       '-o ' + new_output)

    # submit
    q_sub([SAM_commandline], out=new_output.replace('.vcf', ''), t=24*3, jid=job_name, scheduler='SLURM')
