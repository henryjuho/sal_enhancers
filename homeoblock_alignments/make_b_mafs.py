import sys


def main():

    for maf in sys.stdin:

        maf = maf.rstrip()
        out_maf = maf.replace('.maf', '.b.maf')

        with open(out_maf, 'w') as b_maf:

            for line in open(maf):

                new_line = line.rstrip().replace('salmon', 'salmon_b')
                print(new_line, file=b_maf)


if __name__ == '__main__':
    main()