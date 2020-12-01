import argparse
from qsub import q_sub


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-cds_fa', help='Fasta file with CDS sequences in', required=True)
    args = parser.parse_args()

    for i in (0, 2, 3, 4):

        cmd = ('python ~/sal_bal_sel/annotation/degen_to_bed.py '
               '-cds_fa {} -degen {} | '
               'sort -T /scratch/tuyida/bartonhe/tmp/ -k1,1 -k2,2n | '
               'bedtools merge -c 4 -o distinct | '
               'bgzip -c > /scratch/tuyida/bartonhe/sal_ref/salmo_salar_{}fold.bed.gz'
               '').format(args.cds_fa, i, i)

        q_sub([cmd], out='/users/bartonhe/sal_bal_sel/annotation/{}fold_to_bed'.format(i), 
              scheduler='SLURM', rmem=10, mem=10)


if __name__ == '__main__':
    main()
