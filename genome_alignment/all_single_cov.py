from qsub import q_sub
import sys


def main():

    for pairwise_dir in sys.stdin:

        in_dir = pairwise_dir.rstrip()

        single_cov = 'python ~/sal_bal_sel/genome_alignment/single_cov.py -dir {} -ref_name BrownTrout'.format(in_dir)

        q_sub([single_cov], out=in_dir + 'single_cov', scheduler='SLURM')


if __name__ == '__main__':
    main()
