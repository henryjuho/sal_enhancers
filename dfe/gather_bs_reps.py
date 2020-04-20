#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
import sys


def main():

    header = ['run', 'imp', 'exit_code',
              'neu_theta_1', 'neu_gamma_1', 'neu_e_1',
              'sel_theta', 'sel_shape', 'sel_scale', 'sel_e',
              'lnL']

    additional = ['converged', 'bounds', 'rep', 'region']

    head_out = header + additional

    print(*head_out, sep=',')

    res_list = []

    for res in sys.stdin:

        res_stem = res.split('/')[-1]
        bs = int(res_stem.split('.bsrep')[-1].split('.')[0])
        region = res_stem.split('_')[1]

        res = res.rstrip()

        results = an.ResultsFile(open(res))
        mle = results.ml_estimate()
        converged = results.converged()
        bounds = results.bounds_hit(gamma_r=(-500, 100), theta_r=(1e-14, 0.1), r_r=(0.01, 100),
                                    scale_r=(0.1, 5000.0))

        current_res = [mle[x] for x in header] + [converged, ';'.join(bounds), bs, region]

        res_list.append([bs, current_res])

    for processed_res in sorted(res_list, key=lambda x: x[0]):

        print(*processed_res[1], sep=',')


if __name__ == '__main__':
    main()
