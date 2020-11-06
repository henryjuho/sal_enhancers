#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
import argparse
from qsub import q_sub
import os
import sys


def sel_v_neu_anavar(sfs_dat, constraint, n, c, dfe, alg, nnoimp, maximp,
                     out_stem, search, degree, spread, start_index, given):

    """
    submits anavar jobs to cluster after writing required files etc
    :param sfs_dat: dict
    :param constraint: str
    :param n: int
    :param c: int
    :param dfe: str
    :param alg: str
    :param nnoimp: int
    :param maximp: int
    :param out_stem: str
    :param search: int
    :param degree: int
    :param spread: int
    :param start_index: int
    :param given: bool
    :return: None
    """

    anavar_path = ''

    anavar_cmd = '{path}anavar {ctl} {rslts} {log} {seed}'

    # sort file names
    ctl_name = out_stem + '.control.txt'
    merge_out = out_stem + '.merged.results.txt'

    # catch given on first run
    init = ()
    if given:
        if not os.path.isfile(merge_out):
            sys.exit('Given True but no previous runs completed to take besty res from')
        else:
            # get best result from merged out
            best_res = an.ResultsFile(open(merge_out)).ml_estimate(as_string=True)
            init = tuple(best_res.split()[3:-1])

    # make control file
    ctl = an.SNPNeuSelControlFile()

    ctl.set_alg_opts(search=search, alg=alg, key=3,
                     epsabs=1e-20, epsrel=1e-9, rftol=1e-9,
                     maxtime=3600, optional=True,
                     maximp=maximp, nnoimp=nnoimp, init=init)

    ctl.set_data(sfs_dat, n, dfe=dfe, c=c, gamma_r=(-500, 100), theta_r=(1e-14, 0.1), r_r=(0.01, 100),
                 scale_r=(0.001, 10000.0), shape_r=(1e-5, 200), snp_fold=False)
    if degree != 50:
        ctl.set_dfe_optional_opts(degree=degree, optional=True)
    ctl.set_constraint(constraint)
    ctl_contents = ctl.construct()
    with open(ctl_name, 'w') as control:
        control.write(ctl_contents)

    res_file_list = out_stem + '.allres.list.txt'
    with open(res_file_list, 'a') as res_list:

        # split into requested jobs
        for i in range(start_index, start_index+spread):

            split_stem = '{}.split{}'.format(out_stem, i)

            result_name = split_stem + '.results.txt'
            log_name = split_stem + '.log.txt'

            print(result_name, file=res_list)

            # call anavar
            rep_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name, seed=i)

            q_sub([rep_cmd], out=split_stem, jid=split_stem.split('/')[-1] + '.sh', t=72, scheduler='SLURM')


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help='Sample size', required=True, type=int)
    parser.add_argument('-c', help='Number of classes to run model with', required=True, type=int)
    parser.add_argument('-dfe', help='type of dfe to fit, discrete or continuous', default='discrete',
                        choices=['discrete', 'continuous'])
    parser.add_argument('-constraint', help='Constraint for model', choices=['none', 'equal_mutation_rate'],
                        default='none')
    parser.add_argument('-n_search', help='Number of searches to conduct per job', default=500, type=int)
    parser.add_argument('-alg', help='Algorithm to use', default='NLOPT_LD_SLSQP',
                        choices=['NLOPT_LD_SLSQP', 'NLOPT_LD_LBFGS', 'NLOPT_LN_NELDERMEAD',
                                 'NLOPT_LD_TNEWTON_PRECOND_RESTART'])
    parser.add_argument('-nnoimp', help='nnoimp value', default=1, type=int)
    parser.add_argument('-maximp', help='maximp value', default=3, type=int)
    parser.add_argument('-split', help='Number of jobs to split runs across, each job will run the control file once'
                                       'with a different seed given to anavar', default=1, type=int)
    parser.add_argument('-degree', help='changes degree setting in anavar', default=50, type=int)
    parser.add_argument('-out_pre', help='File path and prefix for output', required=True)
    parser.add_argument('-start_index', help='ID for first bin and for first seed', default=1, type=int)
    parser.add_argument('-given', help='If specified takes prev best result as starting values for all runs',
                        default=False, action='store_true')
    parser.add_argument('-target', help='selected sites to use', default='0fold')
    args = parser.parse_args()

    for line in sys.stdin:

        if line.startswith('bs'):
            continue

        # get sfs data from line in right format
        bs, neu_sfs, neu_call, sel_sfs, sel_call = line.rstrip().split('\t')
        neu_sfs = [int(x) for x in neu_sfs.split(',')]
        sel_sfs = [int(x) for x in sel_sfs.split(',')]
        neu_m, sel_m = int(neu_call), int(sel_call)

        sfs_data = {'selected_SNP': (sel_sfs, sel_m), 'neutral_SNP': (neu_sfs, neu_m)}

        out = args.out_pre + '.bsrep' + bs

        # construct process
        sel_v_neu_anavar(sfs_dat=sfs_data,
                         constraint=args.constraint,
                         n=args.n, c=args.c, dfe=args.dfe,
                         alg=args.alg,
                         nnoimp=args.nnoimp, maximp=args.maximp,
                         out_stem=out,
                         search=args.n_search,
                         degree=args.degree,
                         spread=args.split,
                         start_index=args.start_index,
                         given=args.given)


if __name__ == '__main__':
    main()