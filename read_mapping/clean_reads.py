import argparse
import os
import subprocess
from qsub import q_sub


def find_complete(clean_dir):

    clean_files = [x for x in os.listdir(clean_dir) if x.endswith('.fq.gz')]
    val2_samples = set([x.split('_')[0] for x in clean_files if x.endswith('val_2.fq.gz')])
    incomplete_samples = set([x.split('_')[0] for x in clean_files if x.endswith('trimmed.fq.gz')])
    completed_samples = val2_samples - incomplete_samples

    return completed_samples


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-fastq_dir', help='directroy containg fastq files', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    args = parser.parse_args()

    # check out dir for complete cleaned samples
    complete = find_complete(args.out_dir)

    # pair reads
    reads = {}
    for file_name in os.listdir(args.fastq_dir):
        if not file_name.endswith('.fastq.gz'):
            continue

        sample = file_name.split('_')[0]

        if sample in complete:
            continue

        if sample not in reads.keys():
            reads[sample] = []

        reads[sample].append(args.fastq_dir + file_name)

    # cleaning
    for s in reads.keys():

        r1, r2 = sorted(reads[s])

        cmd = ('trim_galore --fastqc --output_dir {out} --paired {r1} {r2}'
               '').format(out=args.out_dir, r1=r1, r2=r2)

        print('running: ' + cmd)

        q_sub([cmd], out=args.out_dir + s, rmem=8, mem=8, scheduler='SLURM', t=2)


if __name__ == '__main__':
    main()
