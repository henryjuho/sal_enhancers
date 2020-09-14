import sys


def main():

    coverage = {i: 0 for i in range(0, 9)}
    counter = 0

    for line in sys.stdin:

        line = line.split('\t')
        bases = int(line[2]) - int(line[1])
        n_cov = 8 - line[5].count('?')
        counter += bases
        coverage[n_cov] += bases

    print('n_species_aligned', 'n_bases', 'percent_alignment', sep=',')
    print('salmon_ref', counter, 1.0, sep=',')

    for i in range(0, 9):

        print(i, coverage[i], round(coverage[i]/counter, 3), sep=',')


if __name__ == '__main__':
    main()
