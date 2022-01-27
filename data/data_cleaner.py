"""
This script provides functions which process GISAID spike protein fasta datafiles to produce
cleaner, more manageable datasets.
"""
from Bio import SeqIO
from collections import defaultdict
import os

script_dir = os.path.dirname(os.path.realpath(__file__))  # path to this file
data_dir = script_dir + '/spike_proteins/'  # relative path of datasets


def reduce_fasta_to_unique_sequences_fasta(infilename, outfilename, data_directory=data_dir):
    """
    Reads a fasta file, identifies the unique sequences and outputs these to new fasta file.
    Also constructs dictionaries of the counts of each unique sequence (sequence_count_data)
    as well as the of number of sequences of each length (sequence_length_dict).
    The new fasta file contains in its description a count of how common the sequence was.


    :param infilename: The input fasta file to be processed.
    :param outfilename: The name of the outputfile.
    :param data_directory: Relative path to the data directory.
    :return: Dictionary counting the number of each type of sequence and the number of
             of each length of sequence.
    """
    fasta_sequences = SeqIO.parse(open(data_directory + infilename), 'fasta')
    sequence_count_dict = defaultdict(lambda: 0)
    sequence_length_dict = defaultdict(lambda: 0)
    number_of_sequences = 0
    for fasta in fasta_sequences:
        identifier, sequence = fasta.id, str(fasta.seq)[:-1]  # remove the last character from the sequence (which is an '*', not an amino acid)
        sequence_count_dict[sequence] += 1
        sequence_length_dict[len(sequence)] += 1
        number_of_sequences += 1

    # sort sequences in decreasing order of frequency of occurrence
    sequence_count_dict = dict(sorted(sequence_count_dict.items(), key=lambda item: item[1], reverse=True))
    with open(data_dir + outfilename, "w") as f:
        for seq,count in sequence_count_dict.items():
            print(f'>{count}', file=f)
            print(seq, file=f)

    print(f'Processed {infilename}:')
    print(f'Found {number_of_sequences} sequences,')
    print(f'of which {len(sequence_count_dict)} are unique sequences.')

    return sequence_count_dict, sequence_length_dict


def remove_corrupt_sequences_from_fasta(infilename, outfilename, length_cutoff=1200, invalid_amino_acids_cutoff=1):
    """
    Given an input fasta file, creates a new fasta file with corrupt sequences removed.
    Corrupt sequences are those which are either too short, or contain too many 'X's, indicating low data quality.

    :param infilename: The file to be read.
    :param outfilename:  The file to be created.
    :param length_cutoff:  Sequences below this length are removed.
    :param invalid_amino_acids_cutoff:  Sequences with this many or more 'X' amino acids are removed.
    """

    with open(data_dir + outfilename, "w") as outfile:
        with open(data_dir + infilename, "r") as infile:
            for sequence in infile:
                # if the line is not a descriptor line
                if sequence[0] != '>':
                    # if sequence is not corrupt
                    if len(sequence) > length_cutoff and sequence.count('X') < invalid_amino_acids_cutoff:
                        outfile.write(descriptor_line)
                        outfile.write(sequence)
                    else:
                        pass
                # temporarily store the descriptor
                else:
                    descriptor_line = sequence


if __name__ == "__main__":
    sequence_counts, sequence_lengths = reduce_fasta_to_unique_sequences_fasta(infilename='1_in_500_spikeprot0112.fasta',
                                                                               outfilename='1_in_500_unique_sequences.fasta')
    remove_corrupt_sequences_from_fasta(infilename='1_in_500_unique_sequences.fasta', outfilename='1_in_500_cleaned.fasta')
