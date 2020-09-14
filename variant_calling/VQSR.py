#!/usr/bin/env python

import argparse
from qsub import q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='VCF file for VQSR to be conducted on', required=True)
    parser.add_argument('-ref', help='Reference genome location', required=True)
    parser.add_argument('-truth_set', help='Truth set of known variants for VQSR', required=True)
    parser.add_argument('-mode', help='Mode to run in', choices=['SNP', 'INDEL'], required=True)
    parser.add_argument('-out', help='Output_directory', required=True)
    args = parser.parse_args()

    # variables
    target_vcf = args.vcf
    ref_genome = args.ref
    truth_set = args.truth_set
    output_prefix = args.out + target_vcf[target_vcf.rfind('/')+1:].rstrip('.vcf')
    recal_file = output_prefix + '.recal'
    tranche_file = output_prefix + '.tranches'
    tranche_list = ['100.0', '99.9', '99.5', '99.0', '98.0', '95.0', '90.0']
    mode = args.mode

    # variant recalibration
    recalibration_commandline = ('gatk VariantRecalibrator '
                                 '-R ' + ref_genome + ' '
                                 '-V ' + target_vcf + ' '
                                 '-resource:hardfilters,known=true,training=true,truth=true,prior=12.0 ' + truth_set + ' '
                                 '-an DP '
                                 '-an QD '
                                 '-an MQ '
                                 '-an MQRankSum '
                                 '-an ReadPosRankSum '
                                 '-an FS '
                                 '-an SOR '
                                 '-an InbreedingCoeff '
                                 '-mode ' + mode + ' '
                                 '-O ' + recal_file + ' '
                                 '--tranches-file ' + tranche_file + ' '
                                 '--rscript-file ' + output_prefix + '.plots.R ')
    for tranche in tranche_list:
        recalibration_commandline += '-tranche ' + tranche + ' '

    # apply recalibration and extract passed
    commandlines_list = [recalibration_commandline, 'wait']
    apply_commandlines_list = []
    extract_commandlines_list = []

    for tranche in tranche_list:
        tranche_out_prefix = output_prefix + '.recalibrated.filtered_t' + tranche

        apply_commandline = ('gatk ApplyVQSR '
                             '-R ' + ref_genome + ' '
                             '-V ' + target_vcf + ' '
                             '--truth-sensitivity-filter-level ' + tranche + ' '
                             '--tranches-file ' + tranche_file + ' '
                             '--recal-file ' + recal_file + ' '
                             '-mode ' + mode + ' '
                             '-O ' + tranche_out_prefix + '.vcf &')

        extract_commandline = ('gatk SelectVariants '
                               '-R ' + ref_genome + ' '
                               '-V ' + tranche_out_prefix + '.vcf '
                               '-O ' + tranche_out_prefix + '.pass.vcf '
                               '--exclude-filtered '
                               '--select-type-to-include ' + mode + ' &')

        apply_commandlines_list.append(apply_commandline)
        extract_commandlines_list.append(extract_commandline)

    commandlines_list = commandlines_list + apply_commandlines_list + ['wait'] + extract_commandlines_list + ['wait']

    # submit jobs
    q_sub(commandlines_list, out=output_prefix + '.recal',
          tr=len(tranche_list), mem=15, rmem=15, scheduler='SLURM')


if __name__ == '__main__':
    main()
