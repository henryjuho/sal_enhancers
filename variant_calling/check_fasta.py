import pysam
import sys


def main():

    fasta = sys.argv[1]
    fasta_true = pysam.FastaFile(sys.argv[2])

    fa_index = fasta + '.fai'

    fasta = pysam.FastaFile(fasta)

    for line in open(fa_index):

        if line.startswith('NC_02'):

            chromo = line.split()[0]

            length = len(fasta.fetch(chromo))
            length2 = len(fasta_true.fetch(chromo))
            print(chromo, length, length2, length == length2, sep='\t')


if __name__ == '__main__':
    main()
