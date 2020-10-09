import argparse
import pysam
from qsub import q_sub


def angsd_region_file(pysam_bed, chromo, out):

    with open(out, 'w') as reg_file:
        for line in pysam_bed.fetch(chromo, parser=pysam.asTuple()):

            new_line = '{}:{}-{}'.format(*line[0:3])
            print(new_line, file=reg_file)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-bam_list', help='list of bam files', required=True)
    parser.add_argument('-ref', help='reference genome to pass to -anc', required=True)
    parser.add_argument('-autosome_list', help='List of autosomes', required=True)
    parser.add_argument('-bed', help='bed file of regions to get sfs for', required=True)
    parser.add_argument('-out', help='Ouput directory and file stem', required=True)
    args = parser.parse_args()

    bed = pysam.TabixFile(args.bed)

    for chromo in open(args.autosome_list):

        chromo = chromo.rstrip()
        chromo_out = args.out + '_' + chromo
        reg_file_name = chromo_out + '.regions.txt'
        angsd_region_file(bed, chromo, reg_file_name)

        angsd_cmd = ('angsd -bam {} '
                     '-doSaf 1 -fold 1 '
                     '-anc {} -ref {} '
                     '-GL 1 -baq 1 -minMapQ 20 -minQ 25 -minInd 31 -setMinDepth 165 -setMaxDepth 372 '
                     '-rf {} -doCounts 1 '
                     '-out {}').format(args.bam_list, args.ref, args.ref, reg_file_name, chromo_out)

        sfs_cmd = 'realSFS {stem}.saf.idx -maxIter 100 > {stem}.sfs'.format(stem=chromo_out)

        q_sub([angsd_cmd, sfs_cmd], out=chromo_out, t=16, rmem=30, mem=30, scheduler='SLURM')


if __name__ == '__main__':
    main()
