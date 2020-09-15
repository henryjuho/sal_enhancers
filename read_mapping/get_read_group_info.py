import sys
import time
import gzip


def sample_names(sample_csv, download_info):

    accession2sample = {x.split(',')[-1].replace(' ', '').rstrip(): x.split(',')[3].replace(' ', '').rstrip() for x in open(sample_csv) if not x.startswith('Pop')}

    run_accession2sample = {}

    for line in open(download_info):

        if line.startswith('study'):
            continue

        line = line.rstrip().split('\t')
        sample_accession = line[2]
        run_accession = line[4]
        try:
            sample = accession2sample[sample_accession]
        except KeyError:
            sample = accession2sample[run_accession]

        run_accession2sample[run_accession] = sample

    return run_accession2sample


def main():

    """
    from Padraic's: https://github.com/padraicc/Corcoran_et_al_2017/blob/master/great_tit/get_read_group_info_bgi.py
    read group tag info: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    also here: https://samtools.github.io/hts-specs/SAMv1.pdf
    """

    # generate time for dt read group tag
    # should be date of sequencing, but just use current
    # date to fill in this tag, so that the file produced
    # works with the bwa_mem_real.sh script for mapping
    date = time.localtime()
    dt = str(date[0]) + '-' + str(date[1]) + '-' + str(date[2])

    get_sample = sample_names('barson_tableS5.csv', 'PRJEB10744.txt')

    f1 = open('fastq_and_rg.txt', 'w')
    for i in sys.stdin:

        i = i.rstrip()
        col = i.strip().split('/')
        run_accession = col[-1].split('_')[0]

        sm = get_sample[run_accession]  # sample acquired from csv file

        # get the flowcell and lane id from the first record of fastq 1
        fq_1 = i
        fq_2 = i.replace('1_val_1.fq.gz', '2_val_2.fq.gz')

        fq = gzip.open(fq_1, 'rb')
        fq_header = fq.readline().decode('utf-8')
        fq.close()

        fq_fields = fq_header.split(':')

        flowcell = fq_fields[2]
        lane = fq_fields[3]

        pl = 'illumina'  # platform
        pu = flowcell + '_' + lane  # platform unit
        lb = run_accession  + '_' + 'lib1'  # library identifier
        pi = '500'  # Predicted median insert size
        ds = 'paired_end'  # description
        cn = 'NORWEGIAN_UNIVERSITY_OF_LIFE_SCIENCES'  # Name of sequencing center producing the read.

        rg_id = lb + '_' + lane  # read group id

        # print fq_1 + '\t' + fq_2 + '\t' + rg_id + '\t' + pl \
        # + '\t' + pu + '\t' + lb + '\t' + pi + '\t' + \
        # ds + '\t' + dt + '\t' + sm + '\t' + cn

        f1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                 (fq_1, fq_2, rg_id, pl, pu, lb, pi, ds, dt, sm, cn))


if __name__ == '__main__':
    main()
