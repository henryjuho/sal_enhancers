import sys
from get_read_group_info import sample_names


get_sample = sample_names('barson_tableS5.csv', 'PRJEB10744.txt')
get_sample = {get_sample[x]: x for x in get_sample.keys()}

out = '| {} | {} | {} | {} |\n|:--:|:--:|:--:|:--:|'.format('run accession', 'sample', 'mean coverage', 'median coverage')
print(out)

for cov_report in sys.stdin:

    cov_report = cov_report.rstrip()

    run_accession = cov_report.split('/')[-1].split('.')[0]
    sm = get_sample[run_accession]

    cov = open(cov_report).readlines()[7].split('\t')
    mean_cov = cov[1]
    med_cov = cov[3]

    out = '| {} | {} | {} | {} |'.format(sm, run_accession, mean_cov, med_cov)
    print(out)

