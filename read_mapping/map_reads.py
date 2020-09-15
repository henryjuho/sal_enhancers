import argparse
import os
from qsub import q_sub, q_write


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-read_info_file', help='text file of fastq pairs and read group info', required=True)
    parser.add_argument('-ref', help='reference file', required=True)
    parser.add_argument('-bam_dir', help='output directory to write bams to', required=True)
    parser.add_argument('-no_sub', help='If specified doesnt submit jobs', default=False, action='store_true')
    args = parser.parse_args()

    for line in open(args.read_info_file):

        fq1, fq2, rg_id, rg_pl, rg_pu, rg_lb, rg_pi, rg_ds, rg_dt, rg_sm, rg_cn = line.rstrip().split('\t')

        bam_id = args.bam_dir + fq1.split('/')[-1].split('_')[0]

        tmp_dir = bam_id + '_tmp'

        # create the temp directory needed by some of the picard tool steps
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)

        bwa_command = ("bwa mem -t 10 -M {ref} {f1} {f2} | "
                       "samtools view -b - > {bam_stem}.bam"
                       "").format(ref=args.ref, f1=fq1, f2=fq2, bam_stem=bam_id)

        # Adds the read group info
        add_read_group = ("java -Xmx12g -jar ~/picard.jar AddOrReplaceReadGroups "
                          "I={bam_stem}.bam O={bam_stem}.id.bam "
                          "RGID={id} RGPL={pl} RGPU={pu} "
                          "RGLB={lb} RGPI={pi} RGDS={ds} "
                          "RGDT={dt} RGSM={sm} RGCN={cn} "
                          "SORT_ORDER=coordinate TMP_DIR={tmp_dir} "
                          "MAX_RECORDS_IN_RAM=2000000"
                          "").format(bam_stem=bam_id, id=rg_id, pl=rg_pl, pu=rg_pu, lb=rg_lb, pi=rg_pi,
                                     ds=rg_ds, dt=rg_dt, sm=rg_sm, cn=rg_cn, tmp_dir=tmp_dir)

        # sort bam file - should be done bove - but rick has this in for some reason?
        # sort_bam = ("java -Xmx12g -jar picard.jar SortSam "
        #            "I={bam_stem}.id.bam "
        #            "O={bam_stem}.sorted.bam "
        #            "SORT_ORDER=coordinate "
        #            "TMP_DIR={tmp_dir} MAX_RECORDS_IN_RAM=2000000"
        #            "").format(bam_stem=bam_id, tmp_dir=tmp_dir)

        # Mark duplicates reduced from 64 to 12 - should be able to get away with due to 8x vs 44x coverage?
        mark_dups = ("java -Xmx12g -jar ~/picard.jar MarkDuplicates "
                     "I={bam_stem}.id.bam O={bam_stem}.dedup.bam "
                     "M={bam_stem}.markdup.metrics.txt "
                     "CREATE_INDEX=true TMP_DIR={tmp_dir} MAX_RECORDS_IN_RAM=2000000"
                     "").format(bam_stem=bam_id, tmp_dir=tmp_dir)

        # Perform a QC step after the mark duplicate step
        # run the validation step

        wgs_metrics = ("java -Xmx12g -jar ~/picard.jar CollectWgsMetrics "
                       "I={bam_stem}.dedup.bam "
                       "O={bam_stem}.wgsmetrics_file.txt "
                       "R={ref} INCLUDE_BQ_HISTOGRAM=true"
                       "").format(bam_stem=bam_id, ref=args.ref)

        validate_bam = ("java -Xmx12g -jar ~/picard.jar ValidateSamFile "
                        "I={bam_stem}.dedup.bam "
                        "O={bam_stem}.validation.txt"
                        "").format(bam_stem=bam_id)

        if args.no_sub:
            q_write([bwa_command, add_read_group, mark_dups, wgs_metrics, validate_bam],
                    out=bam_id + '_mapping', t=48, tr=10, mem=14, rmem=14, scheduler='SLURM')
        else:
            q_sub([bwa_command, add_read_group, mark_dups, wgs_metrics, validate_bam],
                  out=bam_id + '_mapping', t=48, tr=10, mem=14, rmem=14, scheduler='SLURM')

if __name__ == '__main__':
    main()
