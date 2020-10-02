import sys


def main():

    merged_sfs = {}

    for sfs_in in sys.stdin:

        sfs = [float(x) for x in open(sfs_in.rstrip()).readlines()[0].split()]
        n = len(sfs) * 2

        for i in range(0, len(sfs)):

            freq = i / n
            count = sfs[i]

            if freq not in merged_sfs.keys():
                merged_sfs[freq] = 0

            merged_sfs[freq] += count

    for f in sorted(merged_sfs.keys()):

        print(merged_sfs[f], f, sep='\t')


if __name__ == '__main__':
    main()
