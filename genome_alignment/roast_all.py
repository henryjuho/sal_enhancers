import argparse
import os
from qsub import q_write, q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-top_dir', help='top level directory of chromo grouped mafs', required=True)
    parser.add_argument('-ref', help='Reference species name', required=True)
    parser.add_argument('-out', help='Output directory', required=True)
    parser.add_argument('-no_sub', default=False, action='store_true')
    args = parser.parse_args()

    # number of runs
    n_runs = len([x for x in os.listdir(args.top_dir) if x.startswith('BrownTrout')])

    # loop through queries
    sh_list = []
    for i in range(1, n_runs, 200):

        roast_wrapper = ('python ~/sal_bal_sel/genome_alignment/roast_fish.py -top_dir {} -ref {} '
                         '-chr_tag $SLURM_ARRAY_TASK_ID').format(args.top_dir, args.ref)

        start = i

        if n_runs - i < 200:
            end = n_runs
        else:
            end = i + 199

        # if runs above 1000, ie i == 1001 or higher, switch new flag to add multiple of 1000 on
        if i >= 1001:
            multiplier = ' -n_thou 1'
            start -= 1000
            end -= 1000
        else:
            multiplier = ' -n_thou 0'

        q_write([roast_wrapper + multiplier],
                args.out + 'multiz_start' + str(i),
                t=8,
                rmem=4, mem=4,
                array=[start, end],
                scheduler='SLURM')

        sh_list.append(args.out + 'multiz_start' + str(i) + '_job.sh')

    # submit control script
    control = 'python ~/sal_bal_sel/genome_alignment/pairwise_control.py -sh ' + ' -sh '.join(sh_list)
    if args.no_sub:
        q_write([control], out=args.out + 'all_multiple',
                t=72, rmem=2, mem=2,
                scheduler='SLURM')
    else:
        q_sub([control], out=args.out + 'all_multiple',
              t=72, rmem=2, mem=2,
              scheduler='SLURM')


if __name__ == '__main__':
    main()
