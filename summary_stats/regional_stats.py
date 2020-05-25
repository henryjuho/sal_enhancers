from __future__ import print_function
from __future__ import division
import math
import sys
import argparse


def theta_w(n, segsites):
    # theta = S/a
    alpha = sum(1.0/z for z in range(1, n))
    theta = float(segsites) / alpha
    return theta


def pi(n, allele_frq_list):
    no_seqs_dif = [(1.0 - raf**2 - (1.0-raf)**2) * (n/(n-1.0)) for raf in allele_frq_list]
    seq_pi = sum(no_seqs_dif)
    return seq_pi


def tajimas_d(n, allele_frq_list):
    segsites = float(len(allele_frq_list))
    little_d = pi(n, allele_frq_list) - theta_w(n, segsites)

    a1 = sum(1.0/z for z in range(1, n))
    a2 = sum(1.0/z**2 for z in range(1, n))

    e1 = (1.0 / a1) * (((n + 1.0) / (3.0 * (n - 1.0))) - (1.0 / a1))
    e2 = (1.0 / (a1**2 + a2)) * \
         (((2.0 * (n**2 + n + 3.0)) / ((9.0 * n) * (n - 1.0))) -
          ((n + 2.0) / (n * a1)) +
          (a2 / a1**2))
    # print(e2)
    vd = (e1 * segsites) + ((e2 * segsites) * (segsites - 1.0))
    big_d = little_d / math.sqrt(vd)

    return big_d


def string2freq(sfs_string, n):

    freq_list = []
    sfs_counts = [int(x) for x in sfs_string.split(',')]

    for i in range(0, len(sfs_counts)):

        n_snps = sfs_counts[i]
        freq = round((i+1) / n, 3)

        freqs_to_add = [freq for _i in range(0, n_snps)]
        freq_list += freqs_to_add

    return freq_list


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-sel', help='if specified use selected sfs', default=False, action='store_true')
    args = parser.parse_args()
    n = 60

    print('bs', 'region', 'theta', 'pi', 'tajimas_d', 'n_call', sep=',')

    for sfs_file in sys.stdin:

        for line in open(sfs_file.rstrip()):

            if line.startswith('bs_rep'):
                continue

            bs, neu_sfs, neu_call, sel_sfs, sel_call = line.rstrip().split('\t')

            if args.sel:
                sfs = string2freq(sel_sfs, n)
                call = int(sel_call)
                region = sfs_file.split('/')[-1].split('_dfe_')[0]

            else:
                sfs = string2freq(neu_sfs, n)
                call = int(neu_call)
                region = '4fold'

            theta = theta_w(n, len(sfs)) / call
            nuc_div = pi(n, sfs) / call
            d = tajimas_d(n, sfs)

            print(bs, region, theta, nuc_div, d, call, sep=',')


if __name__ == '__main__':
    main()
