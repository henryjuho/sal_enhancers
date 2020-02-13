#!/usr/bin/env python

import argparse
import subprocess
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-ref_fa', help='reference fasta', required=True)
parser.add_argument('-ref_name', help='name of reference species', required=True)
parser.add_argument('-query_fa', help='query fasta', required=True)
parser.add_argument('-query_name', help='name of query species', required=True)
parser.add_argument('-out', help='Output directory', required=True)
parser.add_argument('-hspthresh', default='2200')
args = parser.parse_args()

# variables
out_dir = args.out
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
ref_genome = args.ref_fa
ref_name = args.ref_name
query_genome = args.query_fa
query_name = args.query_name

# identify target data and create 2bit files
target_size_file = out_dir + ref_name + '.sizes'
if not os.path.isfile(target_size_file):
    subprocess.call('faSize -detailed ' + ref_genome + ' > ' + target_size_file, shell=True)
target_2bit = out_dir + ref_name + '.2bit'
if not os.path.isfile(target_2bit):
    subprocess.call('faToTwoBit ' + ref_genome + ' ' + target_2bit, shell=True)

# identify query data and create 2bit files
query = query_genome
query_size_file = out_dir + query_name + '.sizes'
if not os.path.isfile(query_size_file):
    subprocess.call('faSize -detailed ' + query + ' > ' + query_size_file, shell=True)
query_2bit = out_dir + query_name + '.2bit'
if not os.path.isfile(query_2bit):
    subprocess.call('faToTwoBit ' + query + ' ' + query_2bit, shell=True)

# determine output prefix
out_prefix = out_dir + ref_name + '.' + query_name

# lastz commandline
lastz_output = out_prefix + '.axt'
lastz_cmd = ('lastz ' +
             target_2bit + '[nameparse=darkspace] ' +
             query_2bit + '[nameparse=darkspace] '
             '--format=axt '
             '--step=19 '
             '--hspthresh=' + args.hspthresh + ' '
             '--inner=2000 '
             '--ydrop=3400 '
             '--gappedthresh=10000 '
             '--scores=/users/bartonhe/lastz_scores.txt '
             '--chain '
             '> ' + lastz_output)

# chaining commandlines
chain_out = out_prefix + '.chain'
axtChain_cmd = ('axtChain '
                '-linearGap=loose ' +
                lastz_output + ' ' +
                target_2bit + ' ' +
                query_2bit + ' ' +
                chain_out)

prenet_out = out_prefix + '.prenet.chain'
chainPreNet_cmd = ('chainPreNet ' +
                   chain_out + ' ' +
                   target_size_file + ' ' +
                   query_size_file + ' ' +
                   prenet_out)

# netting commandline
target_net_out = out_prefix + '.target.net'
query_net_out = out_prefix + '.query.net'
chain_net_cmd = ('chainNet ' +
                 prenet_out + ' ' +
                 target_size_file + ' ' +
                 query_size_file + ' ' +
                 target_net_out + ' ' +
                 query_net_out)

syn_target_net_out = out_prefix + '.target.syn.net'
syn_query_net_out = out_prefix + '.query.syn.net'
target_netsyntenic_cmd = ('netSyntenic ' + target_net_out + ' ' + syn_target_net_out)
query_netsyntenic_cmd = ('netSyntenic ' + query_net_out + ' ' + syn_query_net_out)

# converting to maf format for use with multiz
post_net_axt = out_prefix + '.postnet.axt'
netToAxt_cmd = ('netToAxt ' +
                syn_target_net_out + ' ' +
                prenet_out + ' ' +
                target_2bit + ' ' +
                query_2bit + ' ' +
                post_net_axt)

sorting_out = out_prefix + '.postnet.sorted.axt'
axtSort_cmd = ('axtSort ' + post_net_axt + ' ' + sorting_out)

output_maf = out_prefix + '.maf'
axtToMaf_cmd = ('axtToMaf ' +
                sorting_out + ' ' +
                target_size_file + ' ' +
                query_size_file + ' ' +
                output_maf)

# run pipeline
subprocess.call(lastz_cmd, shell=True)
subprocess.call(axtChain_cmd, shell=True)
subprocess.call(chainPreNet_cmd, shell=True)
subprocess.call(chain_net_cmd, shell=True)
subprocess.call(target_netsyntenic_cmd, shell=True)
subprocess.call(query_netsyntenic_cmd, shell=True)
subprocess.call(netToAxt_cmd, shell=True)
subprocess.call(axtSort_cmd, shell=True)
subprocess.call(axtToMaf_cmd, shell=True)
