import argparse
from qsub import q_sub


def sample_reads(info_file):

    sample_read_dict = {}

    for line in open(info_file):
        line = line.rstrip().split('\t')
        sample = line [-2]
        accession = line[0].split('/')[-1].split('_')[0]

        if sample not in sample_read_dict.keys():
            sample_read_dict[sample] = []

        sample_read_dict[sample].append(accession)

    return sample_read_dict
 

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-bam_in', help='input directory of unmerged bams', required=True)
    parser.add_argument('-bam_out', help='output directory', required=True)
    parser.add_argument('-info_file', help='file with read group info', required=True)
    parser.add_argument('-ref', help='reference genome', required=True)
    args = parser.parse_args()
    
    all_read_acc = sample_reads(args.info_file)
    
    for fish in all_read_acc.keys():
        
        accessions = all_read_acc[fish]
        in_bams = ' '.join(['I=' + args.bam_in + x + '.dedup.bam' for x in accessions])
        bam_out = args.bam_out + fish + '.dedup.bam'

        merge_cmd = 'java -Xmx12G -jar ~/picard.jar MergeSamFiles {} O={}'.format(in_bams, bam_out)
        wgs_metrics = ("java -Xmx12g -jar ~/picard.jar CollectWgsMetrics "
                       "I={} "
                       "O={}.wgsmetrics_file.txt "
                       "R={} INCLUDE_BQ_HISTOGRAM=true"
                       "").format(bam_out, args.bam_out + fish,  args.ref)
        #print(merge_cmd)
        #print(wgs_metrics)
        q_sub([merge_cmd, wgs_metrics], out=args.bam_out + fish + '_merge', mem=14, rmem=14, scheduler='SLURM')


if __name__ == '__main__':
    main()

