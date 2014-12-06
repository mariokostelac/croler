#!/usr/bin/env python
import optparse
import random

bases = ('A', 'C', 'G', 'T')
header = '>:114:A016U:1:1:17528:1320 1:N:0:1'

reverse_base = { 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A' }

def main():
    parser = optparse.OptionParser()
    parser.add_option('-e', '--error', default=0.02, help='Error rate (per each nucleotide)')
    parser.add_option('-l', '--length', default=1000,help='Length of genome')
    parser.add_option('-r', '--rlength', default=125, help='Length of each read')
    parser.add_option('-c', '--min_rlength', default=60, help='Minimum read length (reads on ends can be shorter)')
    parser.add_option('-m', '--multiplicity', default=2, help='number of rounds (how much time will program go through genome and cut it')
    parser.add_option('-a', '--min_overlap', default=-50, help='minimum overlap length (if it is negative, it means no overlap can occur)')
    parser.add_option('-b', '--max_overlap', default=50, help='maximum overlap length')

    (options, args) = parser.parse_args()

    min_overlap = int(options.min_overlap)
    max_overlap = int(options.max_overlap)
    read_length = int(options.rlength)
    min_rlength = int(options.min_rlength)
    error = float(options.error)
    genome_length = int(options.length)
    multiplicity = int(options.multiplicity)

    # create genome
    genome = []
    for i in xrange(genome_length):
        genome.append(bases[random.randint(0, 3)])

    # create reversed complement
    reversedc = []
    for i in xrange(genome_length):
        reversedc.append(reverse_base[genome[genome_length - i - 1]])

    # approach to cutting
    for i in xrange(multiplicity):
        # every odd step generate reads from reversed complement
        sequence = genome if i % 2 == 0 else reversedc

        spos = 0 if min_overlap >= 0 else random.randint(0, -min_overlap)
        while spos < genome_length:
            read = sequence[spos:spos + read_length]
            spos += read_length - random.randint(min_overlap, max_overlap)

            if len(read) < min_rlength: continue

            # bring some errors
            for x in xrange(len(read)):
                if random.random() < error:
                    read[x] = bases[random.randint(0, 3)]

            print header
            print ''.join(read)


if __name__ == '__main__':
    main()
