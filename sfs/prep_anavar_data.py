import argparse
from collections import Counter
import random


def sfs2counts(freq_list, n):

    """
    converts list of sfs to a condensed sfs that can be plotted or passed to anavar etc
    :param freq_list: list
    :param n: int
    :return: list
    """

    pos_biallelic_freqs = [round(i / float(n), 3) for i in range(1, int(n))]

    counted_sorted_sfs = sorted(Counter([str(x) for x in freq_list]).most_common(), key=lambda z: z[0])
    sfs_freq_dict = {x[0]: x[1] for x in counted_sorted_sfs}

    counts = []
    for frequency in pos_biallelic_freqs:
        try:
            no_var = sfs_freq_dict[str(frequency)]
        except KeyError:
            no_var = 0  # account for freqs with 0 variants

        counts.append(no_var)

    return counts


def read_sfs_file(sfs_file):

    sfs_dict = {}

    for line in open(sfs_file):

        if line.startswith('region'):
            continue

        region, locus, sfs, n_call = line.rstrip().split('\t')

        sfs = [float(x) for x in sfs.split(',')]
        n_call = int(n_call)

        sfs_dict[locus] = [sfs, n_call]

    return sfs_dict


def sfs_for_locus_list(loci, dict1, dict2):

    neu_sfs = []
    sel_sfs = []

    neu_call = 0
    sel_call = 0

    for locus in loci:

        neu_sfs += [round(x, 3) for x in dict1[locus][0]]
        sel_sfs += [round(x, 3) for x in dict2[locus][0]]

        neu_call += dict1[locus][1]
        sel_call += dict2[locus][1]

    return neu_sfs, neu_call, sel_sfs, sel_call


def resample_replace_loci(locus_list):

    resampled_list = []

    for i in range(0, len(locus_list)):

        locus_selected = locus_list[random.randint(0, len(locus_list)-1)]
        resampled_list.append(locus_selected)

    return resampled_list


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-sfs_ref', help='SFS file', required=True)
    parser.add_argument('-sfs_target', help='SFS file', required=True)
    parser.add_argument('-bs_rep', help='Bootstrap replicates', default=100, type=int)
    parser.add_argument('-n', help='sample size', default=31, type=int)
    args = parser.parse_args()

    sfs_ref = read_sfs_file(args.sfs_ref)
    sfs_target = read_sfs_file(args.sfs_target)

    # get list of loci present for both data sources
    all_loci = list([x for x in sfs_ref.keys() if x in sfs_target.keys()])

    print('bs_rep', 'neu_sfs', 'neu_call', 'sel_sfs', 'sel_call', sep='\t')

    for i in range(0, args.bs_rep+1):

        if i == 0:

            locus_list = all_loci

        else:

            locus_list = resample_replace_loci(all_loci)

        neu_sfs, neu_call, sel_sfs, sel_call = sfs_for_locus_list(locus_list, sfs_ref, sfs_target)

        # convert to anavar form
        neu_counts = ','.join([str(x) for x in sfs2counts(neu_sfs, n=2*args.n)])
        sel_counts = ','.join([str(x) for x in sfs2counts(sel_sfs, n=2*args.n)])

        print(i, neu_counts, neu_call, sel_counts, sel_call, sep='\t')


if __name__ == '__main__':
    main()
