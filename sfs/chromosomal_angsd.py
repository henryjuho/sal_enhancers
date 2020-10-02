import argparse
from qsub import q_sub


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-bam_list', help='list of bam files', required=True)
    parser.add_argument('-ref', help='reference genome to pass to -anc', required=True)
    parser.add_argument('-autosome_list', help='List of autosomes', required=True)
    parser.add_argument('-out', help='Ouput directory and file stem', required=True)
    args = parser.parse_args()

    for chromo in open(args.autosome_list):
        
        chromo = chromo.rstrip()
        chromo_out = args.out + '_' + chromo

        angsd_cmd = ('angsd -bam {} '
                     '-doSaf 1 -fold 1 '
                     '-anc {} '
                     '-GL 2 -minMapQ 1 -minQ 20 '
                     '-r {} '
                     '-out {}').format(args.bam_list, args.ref, chromo, chromo_out)

        sfs_cmd = 'realSFS {stem}.saf.idx -maxIter 100 > {stem}.sfs'.format(stem=chromo_out)

        q_sub([angsd_cmd, sfs_cmd], out=chromo_out, rmem=20, mem=20, scheduler='SLURM')


if __name__ == '__main__':
    main()
