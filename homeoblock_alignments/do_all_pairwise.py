import argparse
import subprocess
import os
from qsub import q_sub
import shutil


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-block_info', help='CSV of homeoblock coords', required=True)
    parser.add_argument('-key', help='CSV of chromosome names and corresponding refseq IDs', required=True)
    parser.add_argument('-sal_ref', help='path to salmon reference', required=True)
    parser.add_argument('-sal_query', help='path to salmon query', required=True)
    parser.add_argument('-pike_query', help='path to pike query', required=True)
    parser.add_argument('-out', help='Output directory', required=True)
    args = parser.parse_args()

    keys = {x.split(',')[1].rstrip(): x.split(',')[0] for x in open(args.key)}

    ref_salmon = args.sal_ref
    salmon_query = args.sal_query
    pike = args.pike_query

    for line in open(args.block_info).readlines()[1:]:

        if line.startswith('Block'):
            continue

        # block info
        no, block, ssa, start, end, ssa_q, start_q, end_q, strand = line.split(',')

        # convert to ref_seq
        chromo = keys[ssa]
        chromo_q = keys[ssa_q]

        out_dir = args.out + block + '/'

        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        # make chromosomal fastas for both salmon - homeo block specific
        ref_chromo_fasta = out_dir + 'salmon.fa'
        cmd = 'samtools faidx {} {}:{}-{} > {}'.format(ref_salmon, 'salmon.' + chromo,
                                                       start, end, ref_chromo_fasta)
        subprocess.call(cmd, shell=True)

        query_chromo_fasta = out_dir + 'salmon_b.fa'
        cmd = 'samtools faidx {} {}:{}-{} > {}'.format(salmon_query, 'salmon_b.' + chromo_q,
                                                       start_q, end_q, query_chromo_fasta)
        subprocess.call(cmd, shell=True)

        pike_query = out_dir + 'pike.fa'
        shutil.copy(pike, pike_query)

        # sub align job
        lastz_job1 = ('~/sal_enhancers/homeoblock_alignments/wholegenome_lastz_chain_net.py '
                      '-ref_name salmon -ref_fa {} '
                      '-query_name salmon_b -query_fa {} '
                      '-out {}').format(ref_chromo_fasta, query_chromo_fasta, out_dir)

        lastz_job2 = ('~/sal_enhancers/homeoblock_alignments/wholegenome_lastz_chain_net.py '
                      '-ref_name salmon -ref_fa {} '
                      '-query_name pike -query_fa {} '
                      '-out {}').format(ref_chromo_fasta, pike_query, out_dir)

        q_sub([lastz_job1, lastz_job2], out=out_dir + block, scheduler='SLURM')


if __name__ == '__main__':
    main()
