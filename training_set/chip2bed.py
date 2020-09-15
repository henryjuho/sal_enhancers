import sys


for line in sys.stdin:

    chromo, end = line.rstrip().split()
    start = int(end) - 1

    print(chromo, start, end, sep='\t')
