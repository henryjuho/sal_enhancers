import os
import argparse
import pysam
import ctypes


def find_all(seq, base):

    pos = -1
    while True:
        pos = seq.find(base, pos + 1)
        if pos == -1:
            break
        yield pos


def seq_replace(seq, indices, char='X'):
    
    seq = ctypes.create_unicode_buffer(seq)
    for i in indices:
        seq[i] = char

    return seq.value


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-fa_dir', help='Directory containing alignments and summary csvs', required=True)
    args = parser.parse_args()

    fas = [x for x in os.listdir(args.fa_dir) if x.endswith('fa')]
    print('block', 'n_salmon', 'n_salmon_b', 'n_pike', 'clean_length', 'percent_usable', sep=',')

    for fasta in fas:

        fa = pysam.FastaFile(args.fa_dir + fasta)
        block = fasta.replace('.fa', '')

        sal = fa.fetch('salmon')
        sal_b = fa.fetch('salmon_b')
        pike = fa.fetch('pike')

        raw_len = len(sal) - sal.upper().count('N') - sal.count('-')
        sal_b_aligned = len(sal_b) - sal_b.upper().count('N') - sal_b.count('-')
        pike_aligned = len(pike) - pike.upper().count('N') - pike.count('-')
        
        print(raw_len)

        # get all '-' and 'N' positions
        trim_positions = set()
        for spp_seq in [sal, sal_b, pike]:
            for base in 'N-':
                found_pos = find_all(spp_seq.upper(), base)
                trim_positions |= set(found_pos)
        
        print(len(trim_positions))

        # clean alignment
        trim_positions = sorted(list(trim_positions))
        print(len(trim_positions))
        sal = seq_replace(sal, trim_positions).replace('X', '')
        print(len(sal))
        sal_b = seq_replace(sal_b, trim_positions).replace('X', '')
        pike = seq_replace(pike, trim_positions).replace('X', '')

        # output summary
        print(block, raw_len, sal_b_aligned, pike_aligned, len(sal), round((len(sal)/raw_len) * 100, 3), sep=',')


if __name__ == '__main__':
    main()
