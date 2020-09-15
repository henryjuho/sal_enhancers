import sys
from qsub import q_sub


def main():

    for line in sys.stdin:
        vcf = line.rstrip()

        index = 'gatk IndexFeatureFile -F ' + vcf

        q_sub([index], out=vcf.replace('.g.vcf', '_indexing'), t=1, scheduler='SLURM')


if __name__ == '__main__':
    main()
