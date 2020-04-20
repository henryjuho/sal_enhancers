import sys
import scipy.integrate as integrate
from scipy.stats import gamma
import numpy as np
import pandas as pd


def g_pdf(x, a, scale):
    return gamma.pdf(x, a=a, scale=scale)


def main():

    retain = ['neu_theta_1', 'neu_gamma_1', 'neu_e_1',
              'sel_theta', 'sel_shape', 'sel_scale', 'sel_e',
              'lnL', 'converged', 'bounds', 'rep', 'region']

    new = ['bin', 'proportion']

    header = retain + new

    print(*header, sep=',')

    for res_file in sys.stdin:

        res_file = res_file.rstrip()

        res_data = pd.read_csv(res_file)

        for index, row in res_data.iterrows():

            out_row = [row[x] for x in retain]

            shape = row['sel_shape']
            scale = row['sel_scale']

            # dfe bins for plotting
            bins = [(0, 1), (1, 10), (10, 100), (100, np.inf)]

            for b in bins:

                label = '{} - {}'.format(b[0], b[1])

                # integrate gamma dist probability density function between bin boundaries
                prop = integrate.quad(g_pdf, b[0], b[1], args=(shape, scale))[0]

                out_line = out_row + [label, prop]

                print(*out_line, sep=',')


if __name__ == '__main__':
    main()
