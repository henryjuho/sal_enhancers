import argparse
import sys


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-key', help='csv file of chromo, refseq pairs', required=True)
    args = parser.parse_args()

    keys = {x.split(',')[1].rstrip(): x.split(',')[0] for x in open(args.key)}
    
    for line in sys.stdin:

        if line.startswith('Bloc'):
            continue

        line = line.split(',')

        chr_1, start_1, end_1 = line[2:5]
        chr_2, start_2, end_2 = line[5:8]

        print(keys[chr_1], start_1, end_1, sep='\t')
        print(keys[chr_2], start_2, end_2, sep='\t')


if __name__ == '__main__':
    main()
