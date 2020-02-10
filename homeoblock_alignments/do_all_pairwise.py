import argparse
import subprocess
import os
from qsub import q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-chr_list', help='list of salmon autosomes', required=True)
    parser.add_argument('-ref_fa', help='reference fasta', required=True)
    parser.add_argument('-ref_name', help='name of reference species', required=True)
    parser.add_argument('-query_fa', help='query fasta', required=True)
    parser.add_argument('-query_name', help='name of query species', required=True)
    parser.add_argument('-out', help='Output directory', required=True)
    args = parser.parse_args()

    for line in open(args.chr_list):

        chromo = line.rstrip().split(',')[0]

        out_dir = args.out + '/' + chromo + '/'

        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        query_name = args.query_name + '.' + chromo

        # make chromosomal fastas
        chromo_fasta = out_dir + '/' + query_name + '.fa'
        cmd = 'samtools faidx {} {} > {}'.format(args.query_fa, query_name, chromo_fasta)
        subprocess.call(cmd, shell=True)

        # sub align job
        lastz_wrapper = ('~/sal_enhancers/homeoblock_alignments/wholegenome_lastz_chain_net.py '
                         '-ref_name {} -ref_fa {} '
                         '-query_name {} -query_fa {} '
                         '-out {}').format(args.ref_name, args.ref_fa, query_name, chromo_fasta, out_dir)

        q_sub([lastz_wrapper], out=out_dir + query_name, scheduler='SLURM')


if __name__ == '__main__':
    main()
