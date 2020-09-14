import argparse
import subprocess


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-fa_list', help='Fasta file to add prefix to headers', required=True)
    args = parser.parse_args()

    for line in open(args.fa_list):

        fa, spp = line.rstrip().split()
        fa = args.fa_list[0:args.fa_list.rfind('/')+1] + fa
        rename_cmd = './fasta_add_header_prefix.py -fa {} -pre {}'.format(fa, spp)
        print(rename_cmd)
        subprocess.call(rename_cmd, shell=True)


if __name__ == '__main__':
    main()
