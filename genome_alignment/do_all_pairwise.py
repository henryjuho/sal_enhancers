import argparse
import subprocess
import os
from qsub import q_write, q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-fa_list', help='list of fasta files and spp names', required=True)
    parser.add_argument('-ref_spp', help='Reference species', required=True)
    parser.add_argument('-out', help='Output directory', required=True)
    parser.add_argument('-no_sub', default=False, action='store_true')
    args = parser.parse_args()

    # read in spp
    ref = []
    queries = []

    for line in open(args.fa_list):
        fa, spp = line.rstrip().split()
        fa = args.fa_list[0:args.fa_list.rfind('/')+1] + fa.replace('.fna.gz', '.rename.fa')

        if spp == args.ref_spp:
            ref = [spp, fa]

        else:
            queries.append([spp, fa])

    # count ref contigs
    grep = 'grep -c ">" ' + ref[1]
    n = int(subprocess.Popen(grep, shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[0])

    # loop through queries
    sh_list = []
    for fish in queries:

        out_dir = args.out + ref[0] + '.' + fish[0] + '/'

        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)

        # identify target data and create 2bit files
        target_size_file = out_dir + ref[0] + '.sizes'
        if not os.path.isfile(target_size_file):
            subprocess.call('faSize -detailed ' + ref[1] + ' > ' + target_size_file, shell=True)
        target_2bit = out_dir + ref[0] + '.2bit'
        if not os.path.isfile(target_2bit):
            subprocess.call('faToTwoBit ' + ref[1] + ' ' + target_2bit, shell=True)

        # identify query data and create 2bit files
        query_size_file = out_dir + fish[0] + '.sizes'
        if not os.path.isfile(query_size_file):
            subprocess.call('faSize -detailed ' + fish[1] + ' > ' + query_size_file, shell=True)
        query_2bit = out_dir + fish[0] + '.2bit'
        if not os.path.isfile(query_2bit):
            subprocess.call('faToTwoBit ' + fish[1] + ' ' + query_2bit, shell=True)

        lastz_wrapper = ('~/sal_bal_sel/genome_alignment/wholegenome_lastz_chain_net.py '
                         '-ref_name {} -ref_fa {} '
                         '-query_name {} -query_fa {} '
                         '-out {} -chr_tag $SLURM_ARRAY_TASK_ID').format(ref[0], ref[1], fish[0], fish[1], out_dir)

        # 400 seems to be array limit - array index also cannot exceed 1000
        for i in range(1, n, 200):

            start = i

            if n - i < 200:
                end = n
            else:
                end = i + 199

            # if runs above 1000, ie i == 1001 or higher, switch new flag to add multiple of 1000 on
            if i >= 1001:
                multiplier = ' -n_thou 1'
                start -= 1000
                end -= 1000
            else:
                multiplier = ' -n_thou 0'

            q_write([lastz_wrapper + multiplier],
                    out_dir + ref[0] + '.' + fish[0] + '_start' + str(i),
                    t=8,
                    rmem=4, mem=4,
                    array=[start, end],
                    scheduler='SLURM')

            sh_list.append(out_dir + ref[0] + '.' + fish[0] + '_start' + str(i) + '_job.sh')

    # submit control script
    control = 'python ~/sal_bal_sel/genome_alignment/pairwise_control.py -sh ' + ' -sh '.join(sh_list)
    if args.no_sub:
        q_write([control], out=args.out + 'all_pairwise',
                t=72, rmem=2, mem=2,
                scheduler='SLURM')
    else:
        q_sub([control], out=args.out + 'all_pairwise',
              t=72, rmem=2, mem=2,
              scheduler='SLURM')


if __name__ == '__main__':
    main()
