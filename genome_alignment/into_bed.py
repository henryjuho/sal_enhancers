import argparse
from qsub import q_write, q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-maf', help='MAF file containg alignment, must be compressed', required=True)
    parser.add_argument('-ref_sp', help='Species to use for coordinates in output bed', required=True)
    parser.add_argument('-ref_sizes', help='Sizes file to extract chromosomes for ref species', required=True)
    parser.add_argument('-out', help='Output directory', required=True)
    parser.add_argument('-no_sub', default=False, action='store_true')
    args = parser.parse_args()

    # number of runs
    n_runs = len([x.split()[0] for x in open(args.ref_sizes)])

    # loop through queries
    sh_list = []
    for i in range(1, n_runs, 200):

        wgabed_wrapper = ('python ~/sal_bal_sel/genome_alignment/convert_to_bed.py '
                          '-maf {} -ref_sp {} -ref_sizes {} -out {} '
                          '-chr_tag $SLURM_ARRAY_TASK_ID').format(args.maf, args.ref_sp, args.ref_sizes, args.out)

        # if runs above 1000, ie i == 1001 or higher, switch new flag to add multiple of 1000 on
        n_thou = int(i / 1000)
        start = i % 1000

        multiplier = ' -n_thou {}'.format(n_thou)

        if n_runs - i < 200:
            end = n_runs
        else:
            end = start + 199

        q_write([wgabed_wrapper + multiplier],
                args.out + 'wgabed_start' + str(i),
                t=8,
                rmem=4, mem=4,
                array=[start, end],
                scheduler='SLURM')

        sh_list.append(args.out + 'wgabed_start' + str(i) + '_job.sh')

    # submit control script
    control = 'python ~/sal_bal_sel/genome_alignment/pairwise_control.py -sh ' + ' -sh '.join(sh_list)
    if args.no_sub:
        q_write([control], out=args.out + 'all_wgabed',
                t=72, rmem=2, mem=2,
                scheduler='SLURM')
    else:
        q_sub([control], out=args.out + 'all_wgabed',
              t=72, rmem=2, mem=2,
              scheduler='SLURM')


if __name__ == '__main__':
    main()
