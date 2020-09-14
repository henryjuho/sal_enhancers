import pysam
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    haplotypes = len(list(vcf.header.samples)) * 2
    for line in str(vcf.header).split('\n')[:-1]:
        if 'NW_' in line:
            continue

        print(line)

    # loop through vcf and calc n failed SNPs per individual
    for line in vcf.fetch():

        an = line.info['AN']
        if an != haplotypes:
            continue

        print(str(line).rstrip())


if __name__ == '__main__':
    main()
