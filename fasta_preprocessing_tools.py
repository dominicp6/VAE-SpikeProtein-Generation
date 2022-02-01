"""
This script provides functions which preprocess fasta datafiles by removing incomplete sequences,
collecting repeats and aligning misaligned sequences with MUSCLE.
"""
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from collections import defaultdict
import os
from tqdm import tqdm

script_dir = os.path.dirname(os.path.realpath(__file__))  # path to this file (DO NOT CHANGE)

path_to_muscle_executable = '/home/dominic/miniconda3/pkgs/muscle-3.8.1551-h7d875b9_6/bin/muscle'

# relative path of datasets
data_dir = script_dir + '/data/spike_protein_sequences/'


def remove_incomplete_sequences_from_fasta(infile,
                                           outfile,
                                           length_cutoff=1200,
                                           invalid_amino_acids_cutoff=1,
                                           data_directory = data_dir):
    """
    Given an input fasta file, creates a new fasta file with corrupt sequences removed.
    Corrupt sequences are those which are either too short, or contain too many 'X's, indicating low data quality.

    :param infile: The file to be read.
    :param outfile:  The file to be created.
    :param length_cutoff:  Sequences below this length are removed.
    :param invalid_amino_acids_cutoff:  Sequences with this many or more 'X' amino acids are removed.
    """

    with open(data_directory + outfile, "w") as outfile:
        with open(data_directory + infile, "r") as infile:
            for sequence in tqdm(infile):
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


def downsample_fasta_file(infile,
                          outfile,
                          downsample_factor = 100,
                          data_directory=data_dir):
    """
    Downsamples a fasta file by extracting only every Nth sequence.

    :param infile: The fasta file to downsample.
    :param outfile: The name for the downsampled fasta file to be saved to disk.
    :param downsample_factor: The factor by which to downsample.
    :param data_directory: Relative path to the data directory.
    :return:
    """
    with open(data_directory + infile, "r") as original_fasta_file:
        with open(data_directory + outfile, "w") as downsampled_fasta_file:
            for line_index, line in tqdm(enumerate(original_fasta_file)):
                # factor of 2 because each sequence also has an id row in a fasta file
                print_every = 2 * downsample_factor
                if line_index % print_every == 0 or (line_index - 1) % print_every == 0:
                    print(line, file=downsampled_fasta_file)


def reduce_to_unique_sequences(infile,
                               outfile,
                               data_directory=data_dir):
    """
    Reads a fasta file, identifies the unique sequences and outputs these to new fasta file.
    Also constructs dictionaries of the counts of each unique sequence (sequence_count_data)
    as well as the of number of sequences of each length (sequence_length_dict).
    The new fasta file contains in its description a count of how common the sequence was.


    :param infile: The input fasta file to be processed.
    :param outfile: The name of the outputfile.
    :param data_directory: Relative path to the data directory.
    :return: Dictionary counting the number of each type of sequence and the number of
             of each length of sequence.
    """
    fasta_sequences = SeqIO.parse(open(data_directory + infile), 'fasta')

    # Check is valid fasta file
    if not fasta_sequences:
        raise ValueError(f'{infile} is not a valid .fasta file')

    sequence_count_dict = defaultdict(lambda: 0)    # dict(sequence (str): counts of sequence (int))
    sequence_length_dict = defaultdict(lambda: 0)   # dict(length_of_sqn (int): counts with this length (int))
    number_of_sequences = 0

    for seq_obj in tqdm(fasta_sequences):
        identifier, sequence = seq_obj.id, str(seq_obj.seq)
        # remove non-sequence character suffixes if they exist
        if not sequence[-1].isalpha() and sequence[-1] != '-':
            sequence = sequence[:-1]
        sequence_count_dict[sequence] += 1
        sequence_length_dict[len(sequence)] += 1
        number_of_sequences += 1

    # sort sequences in decreasing order of frequency of occurrence
    sequence_count_dict = dict(sorted(sequence_count_dict.items(), key=lambda item: item[1], reverse=True))
    with open(data_directory + outfile, "w") as f:
        for seq, count in sequence_count_dict.items():
            print(f'>{count}', file=f)
            print(seq, file=f)

    print(f'Processed {infile}:')
    print(f'Found {number_of_sequences} sequences,')
    print(f'of which {len(sequence_count_dict)} are unique sequences.')

    return sequence_count_dict, sequence_length_dict


def reduce_and_align_sequences(infile: str,
                               outfile: str,
                               reduction_factor: int,
                               length_cutoff=1200,
                               invalid_amino_acids_cutoff=1,
                               data_directory = data_dir,
                               path_to_muscle_executable = path_to_muscle_executable):
    """
    Removes incomplete sequences from a fasta database then downsamples, pools identical sequences and aligns them.

    :param infile: The original fasta file to reduce and align.
    :param outfile: The name of the final fasta file.
    :param reduction_factor: Factor by which to downsample.
    :param length_cutoff: Sequences below this length are removed.
    :param invalid_amino_acids_cutoff: Sequences with at least this many invalid amino acids are removed.
    """
    # 1.
    print('Removing incomplete sequences...')
    remove_incomplete_sequences_from_fasta(infile=infile,
                                           outfile=infile + '.cleaned',
                                           length_cutoff=length_cutoff,
                                           invalid_amino_acids_cutoff=invalid_amino_acids_cutoff,
                                           data_directory=data_dir)

    # 2.
    print(f'Downsampling sequences by factor of {reduction_factor}...')
    downsample_fasta_file(infile=infile + '.cleaned',
                          outfile=infile + '.cleaned.downsampled',
                          downsample_factor=reduction_factor,
                          data_directory=data_dir)

    # 3.
    print('Getting unique sequences...')
    reduce_to_unique_sequences(infile=infile + '.cleaned.downsampled',
                               outfile=infile + '.cleaned.downsampled.unique',
                               data_directory=data_dir)

    # 4.
    print('To complete the alignment step run this command in the terminal:')
    muscle_command = MuscleCommandline(path_to_muscle_executable,
                      input=data_directory+infile+'.cleaned.downsampled.unique',
                      out=data_directory+outfile)
    print(muscle_command)


if __name__ == "__main__":
    reduce_and_align_sequences(infile='spikeprot0112.fasta',
                               outfile='aligned_spike_proteins.fasta',
                               reduction_factor=100,
                               length_cutoff=1200,
                               invalid_amino_acids_cutoff=1)