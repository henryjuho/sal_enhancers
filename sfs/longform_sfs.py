import sys


def main():

    print('bs_rep', 'freq', 'count', 'proportion', 'region', sep=',')

    for txt in sys.stdin:

        txt = txt.rstrip()

        region = txt.split('/')[-1].split('_')[0]

        for line in open(txt):

            if line.startswith('bs'):
                continue

            bs, neu_dfs, neu_call, sel_sfs, sel_call = line.rstrip().split('\t')

            sfs = [int(x) for x in sel_sfs.split(',')]
            sfs_prop = [y / sum(sfs) for y in sfs]

            for i in range(0, len(sfs)):

                print(bs, i+1, sfs[i], sfs_prop[i], region, sep=',')


if __name__ == '__main__':
    main()
