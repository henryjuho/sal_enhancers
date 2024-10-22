#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
import sys
import argparse
import re


def aic(n_param, max_lnl):

    """
    calculates AIC as desicribed here: https://en.wikipedia.org/wiki/Akaike_information_criterion
    :param n_param: int
    :param max_lnl: float
    :return:
    """

    k = n_param
    lnl = max_lnl

    return (2 * k) - (2 * lnl)


def delta_aic(current_d, best_d):

    """
    returns delta_AIC, which is defined as the difference between the AIC for the best fitting model
    and the model on the row of the spreadsheet - from Kai's email 14th Sept 17
    :param current_d: float
    :param best_d: float
    :return:
    """

    return best_d - current_d


def reformat_mle(line, n_classes, converged, b_hit, model, p, dfe, n_search,
                 custom_col=None, custom_val=None):

    """
    converts mle estimate line from anavar results file to multiple lines in long form
    :param line: dict
    :param n_classes: int
    :param converged: bool
    :param b_hit: list
    :param model: str
    :param p: int
    :param dfe: str
    :param n_search: int
    :param custom_col: str
    :param custom_val: str
    :return: list
    """

    # prepare header
    out_lines = []
    new_header = ['run', 'imp', 'exit_code',
                  'theta', 'scale', 'shape', 'gamma', 'e',
                  'site_class', 'sel_type', 'lnL',
                  'boundaries_hit', 'searches', 'converged', 'model', 'params']

    non_variable_vals = {'run': line['run'], 'imp': line['imp'], 'exit_code': line['exit_code'],
                         'lnL': line['lnL'], 'boundaries_hit': '|'.join(b_hit), 'converged': converged,
                         'model': model, 'params': p, 'AIC': aic(p, line['lnL']), 'searches': n_search}

    if custom_col is not None:
        new_header.append(custom_col)
        non_variable_vals[custom_col] = custom_val

    new_header.append('AIC')

    # process neutral variants then selected variants
    for sel in ['neu_', 'sel_']:

        # only ever one class of neutral variants
        if sel == 'neu_':
            n_class = 1
        # if sel then number of classes given
        else:
            n_class = n_classes

        for site_class in range(1, n_class+1):

            # get value for each column in output table
            row = []
            for col in new_header:

                # these values same for each row
                if col in non_variable_vals.keys():
                    col_val = non_variable_vals[col]

                elif col == 'site_class':
                    col_val = site_class

                elif col == 'sel_type':
                    col_val = sel.rstrip('_')

                # leaving 'theta', 'scale', 'shape', 'gamma', 'e'
                else:
                    if dfe == 'continuous' and sel == 'sel_':
                        class_str = ''
                    else:
                        class_str = '_' + str(site_class)

                    key_str = '{}{}{}'.format(sel, col, class_str)

                    try:
                        col_val = line[key_str]
                    except KeyError:
                        col_val = 'na'

                row.append(col_val)

            out_lines.append(row)

    return new_header, out_lines


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-file_pattern', help='takes a regular expression in order to extract a custom ID from a file'
                                              'name, along with a column name eg) degree,_degree(\d+)\.')
    parser.add_argument('-dn', help='dn estimates', type=float, required=True)
    parser.add_argument('-ds', help='ds estimates', type=float, required=True)
    args = parser.parse_args()
    counter = 0

    if args.file_pattern is not None:
        spec_col = args.file_pattern.split(',')[0]
        spec_pattern = args.file_pattern.split(',')[1]
    else:
        spec_col, spec_pattern = None, None

    all_res = []
    for res in sys.stdin:

        # sort custom col contents
        if spec_col is not None:
            spec_val = re.search(spec_pattern, res).group(1)
        else:
            spec_val = None

        res = res.rstrip()

        if 'equal_t' in res:
            constraint = '_equal_t'
        else:
            constraint = ''

        results = an.ResultsFile(open(res))
        alpha = results.get_alpha(dn=args.dn, ds=args.ds, var_type='snp')
        mle = results.ml_estimate()
        n_class = results.num_class()
        variant = results.data_type()
        dfe = results.dfe()
        free_params = len(results.free_parameters())

        mod_name = '{}_{}_{}class{}'.format(variant, dfe, n_class, constraint)

        reformed = reformat_mle(mle, n_class, results.converged(),
                                results.bounds_hit(gamma_r=(-500, 100), theta_r=(1e-14, 0.1), r_r=(0.01, 100),
                                scale_r=(0.1, 5000.0)), mod_name, free_params, dfe, results.num_runs(),
                                custom_col=spec_col, custom_val=spec_val)
        reformed[0].append('alpha')
        if counter == 0:
            all_res.append(reformed[0])
            counter += 1

        for x in reformed[1]:
            x.append(alpha)
            all_res.append(x)

    # calc AIC
    aics = sorted([z[-2] for z in all_res[1:]])
    best_aic = aics[0]

    print(*all_res[0] + ['delta_AIC'], sep=',')
    
    out_data = []
    for line in all_res[1:]:
        delta = delta_aic(line[-2], best_aic)
        out_data.append((delta, line + [delta]))
    
    for line in sorted(out_data, reverse=True, key=lambda x: x[0]):
        print(*line[1], sep=',')


if __name__ == '__main__':
    main()
