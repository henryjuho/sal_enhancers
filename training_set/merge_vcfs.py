import sys
import argparse
from qsub import q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-out', help='Output vcf', required=True)
    args = parser.parse_args()

    vcfs = ' I='.join([x.rstrip() for x in sys.stdin])

    cmd = ('java -Xmx10G -jar /users/bartonhe/picard.jar GatherVcfs '
           'I={} O={}').format(vcfs, args.out)

    q_sub([cmd], out=args.out.replace('.vcf', ''), mem=12, rmem=12, scheduler='SLURM')


if __name__ == '__main__':
    main()
