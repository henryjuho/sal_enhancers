import argparse
import sys
import numpy as np
import scipy.integrate as integrate
from scipy.stats import gamma


def gfix_pdf(x, a, scale):
    return (-x / (1-np.exp(x))) * gamma.pdf(x, a=a, scale=scale)


def alpha_continuous(shape, scale, dn, ds):

    u = integrate.quad(gfix_pdf, 0, np.inf, args=(shape, scale))[0]

    alpha = 1 - ((ds/dn) * u)

    return alpha


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-div', help='divergence estimates', required=True)
    args = parser.parse_args()

    div_data = {x.split('\t')[0].replace('_', '-'): float(x.split('\t')[1])
                for x in open(args.div) if not x.startswith('reg')}
    keys = []

    print('region', 'bs_rep', 'theta', 'shape', 'scale', 'mean_gamma', 'alpha', sep=',')

    for line in sys.stdin:
        line = line.rstrip().split(',')

        if line[0] == 'run':
            keys = line
            continue

        bs_dat = {keys[i]: line[i] for i in range(0, len(line))}

        theta = float(bs_dat['sel_theta'])
        shape = float(bs_dat['sel_shape'])
        scale = float(bs_dat['sel_scale'])
        mean_gamma = shape * scale
        median_gamma = gamma.median(a=shape, scale=scale)

        rep = bs_dat['rep']
        region = bs_dat['region']

        print(theta, mean_gamma, median_gamma, div_data['4fold'], div_data[region])

        line_alpha = alpha_continuous(shape=shape, scale=scale, ds=div_data['4fold'], dn=div_data[region])

        print(region, rep, theta, shape, scale, mean_gamma, line_alpha, sep=',')


if __name__ == '__main__':
    main()




