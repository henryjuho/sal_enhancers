import argparse
from qsub import q_sub


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-bed_regs', help='bed file with subregions to get', required=True)
    parser.add_argument('-bed_targets', help='bed file of target sites', required=True)
    parser.add_argument('-vcf', help='VCF file with called SNPs', required=True)
    parser.add_argument('-call_fa', help='Fasta file of callable sites', required=True)
    parser.add_argument('-region', help='Label to put in region column of output', required=True)
    parser.add_argument('-chromo_list', help='list of chromosomes to process', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    args = parser.parse_args()

    for line in open(args.chromo_list):

        chromo = line.rstrip()

        output_file = '{}sfs_ncall_regional.{}.{}.txt'.format(args.out_dir, args.region, chromo)
        log_file = '{}sfs_ncall_regional.{}.{}.log.txt'.format(args.out_dir, args.region, chromo)

        if args.bed_regs == 'NA':
            cmd = ('python ~/sal_enhancers/dfe/bootstrappable_sfs_data_enhancers.py '
                   '-bed_target {bed_target} '
                   '-vcf {vcf} '
                   '-call_fa {call_fa} '
                   '-region {region} '
                   '-chr {chromo} '
                   '1> {out} '
                   '2> {log}').format(bed_regs=args.bed_regs,
                                      bed_target=args.bed_targets,
                                      vcf=args.vcf,
                                      call_fa=args.call_fa,
                                      region=args.region,
                                      chromo=chromo,
                                      out=output_file,
                                      log=log_file)

        else:
            cmd = ('python ~/sal_enhancers/dfe/bootstrappable_sfs_data.py '
                   '-bed_regs {bed_regs} '
                   '-bed_target {bed_target} '
                   '-vcf {vcf} '
                   '-call_fa {call_fa} '
                   '-region {region} '
                   '-chr {chromo} '
                   '1> {out} '
                   '2> {log}').format(bed_regs=args.bed_regs,
                                      bed_target=args.bed_targets,
                                      vcf=args.vcf,
                                      call_fa=args.call_fa,
                                      region=args.region,
                                      chromo=chromo,
                                      out=output_file,
                                      log=log_file)

        q_sub([cmd], out=output_file.replace('.txt', ''), mem=8, rmem=8, t=48, scheduler='SLURM')


if __name__ == '__main__':
    main()
