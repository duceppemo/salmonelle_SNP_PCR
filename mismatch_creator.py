#!/bin/python

__author__ = 'duceppemo'

'''
Take a fasta file containing primer sequences and create all possible primers shuffling
 all the bases except the last one
'''


class MismatchGenerator(object):

    def __init__(self, args):
        import os
        import sys

        # Define variables based on supplied arguments
        self.args = args
        self.input = args.input
        self.output = args.output
        self.mismatch = args.mismatch

        # Checks
        if not os.path.isfile(self.input):
            sys.exit('Supplied input file does not exists.')
        if not isinstance(self.mismatch, (int)):
            sys.exit('Mismatch value must be an integer.')
        if not 1 <= self.mismatch <= 3:
            sys.exit('Mismatch value must be equal to 1, 2 or 3.')

        # create dictionaries to hold data
        self.primer_dict = dict()

        # run the script
        self.run()

    def run(self):
        self.parse_input_fasta()
        self.create_primers()
        self.write_output_fasta()

    def parse_input_fasta(self):
        with open(self.input, 'r') as f:  # open file
            for line in f:  # read file line by line
                line = line.rstrip()  # chomp -> remove trailing whitespace characters
                if line:  # skip blank lines or lines with only whitespaces
                    if line.startswith('>'):  # header
                        [position,  allele] = line.split('_')
                        position = position.strip('>')
                    else:
                        sequence = line
                        self.primer_dict.setdefault(position, {}).setdefault(allele, {}).setdefault(sequence, {})

    def create_primers(self):
        alphabet = ['A', 'T', 'C', 'G']

        for position, alleles in iter(self.primer_dict.items()):
            for allele, sequences in iter(alleles.items()):
                for seq in sequences:
                    snp = seq[-1:]
                    # Create all possible combination of SNPs at each position (1 mismatch)
                    end = len(seq) - 1
                    if self.mismatch == 1:  # 17 primers per primer for a 21-mer
                        for i in range(0, end, 1):  # range([start], stop[, step])
                            for letter1 in alphabet:
                                s1 = list(seq)  # Convert string to list of character
                                s1[i] = letter1  # Replace at desired index
                                m1 = "".join(s1)  # Remake string by joining all element of the list
                            # Write new string to file and append ech variant with a different number
                            self.primer_dict.setdefault(position, {}).setdefault(allele, {}) \
                                .setdefault(seq, {}).setdefault(m1, [])
                    elif self.mismatch == 2:  # 1771 primers per primer for a 21-mer
                        for i in range(0, end, 1):  # range([start], stop[, step])
                            for letter1 in alphabet:
                                s1 = list(seq)  # Convert string to list of character
                                s1[i] = letter1  # Replace at desired index
                                m1 = "".join(s1)  # Remake string by joining all element of the list

                                if i < end:
                                    for j in range((i + 1), end, 1):
                                        for letter2 in alphabet:
                                            s2 = list(m1)
                                            s2[j] = letter2
                                            m2 = "".join(s2)

                                            self.primer_dict.setdefault(position, {}).setdefault(allele, {}) \
                                                .setdefault(seq, {}).setdefault(m2, [])
                    elif self.mismatch == 3:  # 32551 primers per primer for a 21-mer
                        for i in range(0, end, 1):  # range([start], stop[, step])
                            for letter1 in alphabet:
                                s1 = list(seq)  # Convert string to list of character
                                s1[i] = letter1  # Replace at desired index
                                m1 = "".join(s1)  # Remake string by joining all element of the list

                                if i < end:
                                    for j in range((i + 1), end, 1):
                                        for letter2 in alphabet:
                                            s2 = list(m1)
                                            s2[j] = letter2
                                            m2 = "".join(s2)

                                            if j < end:
                                                for k in range((j + 1), end, 1):
                                                    for letter3 in alphabet:
                                                        s3 = list(m2)
                                                        s3[k] = letter3
                                                        m3 = "".join(s3)

                                                        self.primer_dict.setdefault(position, {})\
                                                                        .setdefault(allele, {})\
                                                                        .setdefault(seq, {})\
                                                                        .setdefault(m3, [])
                                                        test = 1

    def write_output_fasta(self):
        fh = open(self.output, 'w')

        for position, alleles in iter(self.primer_dict.items()):
            for allele, sequences in iter(alleles.items()):
                counter = 0
                for sequence, all_primers in iter(sequences.items()):
                    for p in all_primers:
                        counter += 1
                        fh.write('>' + position + '_' + allele + '_' + str(counter) + "\n" + p + "\n")


if __name__ == '__main__':

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Generate all possible degenerated primers from a list of primers')
    parser.add_argument('-i', '--input', metavar='primers.fasta',
                        required=True,
                        help='Input primer fasta file')

    parser.add_argument('-o', '--output', metavar='degenerated_primers.fasta',
                        required=True,
                        help='Output degenerated primer files')

    parser.add_argument('-m', '--mismatch', metavar='1', type=int,
                        required=True,
                        help='Number of mismatched. 1-3')

    # Get the arguments into an object
    arguments = parser.parse_args()

    MismatchGenerator(arguments)
