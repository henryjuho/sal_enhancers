import pysam
import sys
import subprocess


def main():

    fastas = [x.rstrip() for x in sys.stdin]

    for fa in fastas:

        chromo = subprocess.Popen('head -n 1 ' + fa, shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split('\n')[0]
        chromo = chromo.replace('>', '')
        start, stop = chromo.split(':')[1].split('-')
        head_length = int(stop) - int(start)

        fasta = pysam.FastaFile(fa)
        length = len(fasta.fetch(chromo))

        print(fa, chromo, head_length, length, head_length==length, sep='\t')


if __name__ == '__main__':
    main()
