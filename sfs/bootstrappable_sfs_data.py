import argparse
import gzip
import subprocess
import sys


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-bed_regs', help='bed file with subregions to get', required=True)
    parser.add_argument('-bed_target', help='bed file of target sites', required=True)
    parser.add_argument('-vcf', help='VCF file with called SNPs', required=True)
    parser.add_argument('-call_fa', help='Fasta file of callable sites', required=True)
    parser.add_argument('-region', help='Label to put in region column of output', required=True)
    parser.add_argument('-chr', help='chromosome', required=True)
    args = parser.parse_args()

    # loop through bed file
    counter = 0

    print('region', 'window', 'sfs', 'n_call', sep='\t', file=sys.stdout)

    for line in gzip.open(args.bed_regs):
        counter += 1

        # process bed line
        line = line.decode('utf-8').rstrip().split('\t')
        chromo, start, stop = line[:3]

        if chromo != args.chr:
            continue

        if len(line) == 3:
            window = counter
        else:
            window = line[3]

        # get call
        call_cmd = ('tabix {} {}:{}-{} | '
                    'bedtools getfasta -fi {} -bed stdin | '
                    'grep -v ">" | tr -d "\\n"').format(args.bed_target, chromo, start, stop, args.call_fa)

        print(call_cmd, file=sys.stderr)

        n_call = subprocess.Popen(call_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
        n_call = n_call.decode('utf-8').count('K')

        # get sfs
        sfs_cmd = ('tabix -h {} {}:{}-{} | '
                   'bedtools intersect -header -a stdin -b {} | '
                   '~/sfs_utils/vcf2raw_sfs.py '
                   '-mode snp').format(args.vcf, chromo, start, stop, args.bed_target)

        print(sfs_cmd, file=sys.stderr)

        sfs = subprocess.Popen(sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
        sfs = sfs.decode('utf-8').split('\n')[:-1]
        #sfs = [x.split()[0] for x in sfs]
        if len(sfs) == 0:
            sfs = '0'
        else:
            sfs = ','.join(sfs)

        print(args.region, window, sfs, n_call, sep='\t', file=sys.stdout)


if __name__ == '__main__':
    main()
