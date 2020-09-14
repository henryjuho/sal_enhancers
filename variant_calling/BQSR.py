import sys
import argparse
from qsub import q_sub


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-ref', help='reference genome', required=True)
    parser.add_argument('-train_vcf', help='vcf file of training data', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    args = parser.parse_args()

    bams = [x.rstrip() for x in sys.stdin]

    for bam in bams:

        bam_stem = bam.replace('.bam', '').split('/')[-1]
        out_table = '{}{}.table'.format(args.out_dir, bam_stem)
        recal_bam = out_table.replace('.table', '.bqsr.bam')

        bqsr = ('gatk BaseRecalibrator '
                '-I {bam} '
                '-R {ref} '
                '--known-sites {truth} '
                '-O {table}').format(bam=bam, ref=args.ref, truth=args.train_vcf, table=out_table)

        apply = ('gatk ApplyBQSR '
                 '-R {ref} '
                 '-I {bam} '
                 '--bqsr-recal-file {table} '
                 '-O {new_bam}').format(ref=args.ref, bam=bam, table=out_table, new_bam=recal_bam)

        q_sub([bqsr, apply], out=out_table.replace('.table', ''), t=24, scheduler='SLURM')


if __name__ == '__main__':
    main()
