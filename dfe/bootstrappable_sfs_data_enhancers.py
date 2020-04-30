import argparse
import subprocess
import sys


def read_bed_loci(bed_file, chromo):

    bed_cmd = 'zgrep ^{} {} | cut -f 4 | sort -u | grep -v "\\."'.format(chromo, bed_file)

    gs = subprocess.Popen(bed_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[:-1]

    return gs


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-bed_target', help='bed file of target sites', required=True)
    parser.add_argument('-vcf', help='VCF file with called SNPs', required=True)
    parser.add_argument('-call_fa', help='Fasta file of callable sites', required=True)
    parser.add_argument('-region', help='Label to put in region column of output', required=True)
    parser.add_argument('-chr', help='chromosome', required=True)
    args = parser.parse_args()

    # loop through bed file
    counter = 0

    print('region', 'window', 'sfs', 'n_call', sep='\t', file=sys.stdout)

    for locus in read_bed_loci(args.bed_target, args.chr):

        counter += 1

        # get call
        call_cmd = ('zgrep -w {} {} | '
                    'bedtools getfasta -fi {} -bed stdin | '
                    'grep -v ">" | tr -d "\\n"').format(locus, args.bed_target, args.call_fa)

        print(call_cmd, file=sys.stderr)

        n_call = subprocess.Popen(call_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
        n_call = n_call.decode('utf-8').count('K')

        # get sfs
        sfs_cmd = ('zgrep -w {} {} | '
                   'bedtools intersect -header -a {} -b stdin | '
                   '~/sfs_utils/vcf2raw_sfs.py '
                   '-mode snp').format(locus, args.bed_target, args.vcf)

        print(sfs_cmd, file=sys.stderr)

        sfs = subprocess.Popen(sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0]
        sfs = sfs.decode('utf-8').split('\n')[:-1]
        #sfs = [x.split()[0] for x in sfs]
        if len(sfs) == 0:
            sfs = '0'
        else:
            sfs = ','.join(sfs)

        print(args.region, locus, sfs, n_call, sep='\t', file=sys.stdout)


if __name__ == '__main__':
    main()
