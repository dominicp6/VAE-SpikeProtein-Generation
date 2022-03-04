"""
This script provides functions which preprocess fasta datafiles by removing incomplete sequences,
collecting repeats, labelling with variant names, and aligning misaligned sequences with MUSCLE.
"""

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Align.Applications import MuscleCommandline
from collections import defaultdict
import os
import numpy as np
from tqdm import tqdm

script_dir = os.path.dirname(os.path.realpath(__file__))  # path to this file (DO NOT CHANGE)

path_to_muscle_executable = '/home/dominic/miniconda3/pkgs/muscle-3.8.1551-h7d875b9_6/bin/muscle'
path_to_consensus_sequences = os.path.join(script_dir, "data", "spike_protein_sequences", "consensus_sequences")

# relative path of datasets
data_dir = os.path.join(script_dir, "data", "spike_protein_sequences")


def check_all_sequences_have_same_length(database):
    db_seq = SeqIO.parse(database, 'fasta')

    max_length = 0
    error_found = False
    for seq_id, seq in tqdm(enumerate(db_seq)):
        if max_length == 0:
            max_length = len(seq.seq)
        else:
            if len(seq.seq) != max_length:
                print(f"Error: The sequence number {seq_id+1} has length {len(seq.seq)} whereas the first sequence has length {max_length}.")
                error_found = True

    if not error_found:
        print('All sequences have the same length!')


def remove_id_label_from_fasta_database(database, id_position, outfile):
    db_seq = SeqIO.parse(database, 'fasta')

    with open(outfile, 'w') as out_file:
        for seq in tqdm(db_seq):
            seq_id_list = seq.id.split('|')
            del seq_id_list[id_position]
            print(f">{'|'.join(seq_id_list)}", file=out_file)
            print(seq.seq, file=out_file)


def add_id_label_to_fasta_database(database, id_position, label, outfile):
    db_seq = SeqIO.parse(database, 'fasta')

    with open(outfile, 'w') as out_file:
        for seq in tqdm(db_seq):
            seq_id_list = seq.id.split('|')
            seq_id_list.insert(id_position, label)
            print(f">{'|'.join(seq_id_list)}", file=out_file)
            print(seq.seq, file=out_file)


def partition_fasta_database_into_chunks(database, chunk_size, outfile_prefix):
    """
    Splits a single database into multiple smaller databases each containing chunk_size number of sequences
    (except possibly the last chuck, which may be smaller).
    """
    db_seq = SeqIO.parse(database, 'fasta')

    file_index = 0
    saw_a_new_sequence = True
    with tqdm(total=0) as pbar:
        while saw_a_new_sequence:
            for _ in (True,):  # "breakable scope" idiom
                with open(f"{outfile_prefix}_{file_index}.fasta", "w") as out_file:
                    saw_a_new_sequence = False
                    for seq_number, seq in enumerate(db_seq):
                        saw_a_new_sequence = True
                        if seq_number < chunk_size:
                            print(f">{seq.id}", file=out_file)
                            print(seq.seq, file=out_file)
                        else:
                            break
                    file_index += 1
                    pbar.update(1)
                    break


def combine_two_databases(database1, database2, variant_database2, outfile, variant_database1=None):
    """
    Merges two fasta databases, with optional labelling to distinguish sequences originating from the two databases.

    :param database1: Path to first database.
    :param database2: Path to second database.
    :param variant_database2: The "variant" label for the sequences in database 2 (e.g. 'natural', 'VAEsynthetic', ...)
    :param outfile: The name for the outfile of the merged databases.
    :param variant_database1: The "variant" label for the sequences in database 1. If "None", then leave the label
    unchanged.
    """
    db1_seq = SeqIO.parse(database1, 'fasta')
    db2_seq = SeqIO.parse(database2, 'fasta')

    with open(outfile, "w") as out_file:
        for seq in db1_seq:
            if variant_database1 is not None:
                print(f">{seq.id}|{variant_database1}", file=out_file)
            else:
                print(f">{seq.id}", file=out_file)
            print(seq.seq, file=out_file)
        for seq in db2_seq:
            print(f">{seq.id}|{variant_database2}", file=out_file)
            print(seq.seq, file=out_file)


def remove_incomplete_sequences_from_fasta(infile,
                                           outfile,
                                           length_cutoff=1200,
                                           invalid_amino_acids_cutoff=1,
                                           data_directory=data_dir):
    """
    Given an input fasta file, creates a new fasta file with corrupt sequences removed.
    Corrupt sequences are those which are either too short, or contain too many 'X's, indicating low data quality.

    :param infile: The file to be read.
    :param outfile:  The file to be created.
    :param length_cutoff:  Sequences below this length are removed.
    :param invalid_amino_acids_cutoff:  Sequences with this many or more 'X' amino acids are removed.
    :param data_directory: Path to directory of infile and outfile.
    """

    with open(os.path.join(data_directory, outfile), "w") as outfile:
        with open(os.path.join(data_directory, infile), "r") as infile:
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
                          downsample_factor=100,
                          data_directory=data_dir):
    """
    Downsamples a fasta file by extracting only every Nth sequence.

    :param infile: The fasta file to downsample.
    :param outfile: The name for the downsampled fasta file to be saved to disk.
    :param downsample_factor: The factor by which to downsample.
    :param data_directory: Path to directory of infile and outfile.
    """

    with open(os.path.join(data_directory, infile), "r") as original_fasta_file:
        with open(os.path.join(data_directory, outfile), "w") as downsampled_fasta_file:
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
    :param data_directory: Path to directory of infile and outfile.
    :return: Dictionaries counting the number of each type of sequence and the number of
             of each length of sequence.
             e.g. ['ABC' : 1000, 'BCDE' : 500], [3: 1000, 4: 500]
    """
    fasta_sequences = SeqIO.parse(os.path.join(data_directory, infile), 'fasta')

    # Check is valid fasta file
    if not fasta_sequences:
        raise ValueError(f'{infile} is not a valid .fasta file')

    sequence_count_dict = defaultdict(lambda: 0)  # dict(sequence (str): counts of sequence (int))
    sequence_length_dict = defaultdict(lambda: 0)  # dict(length_of_sqn (int): counts with this length (int))
    sequence_date_dict = defaultdict(lambda: [])  # dict(sequence (str): median date that sequence was recorded)
    number_of_sequences = 0

    for seq_obj in tqdm(fasta_sequences):
        identifier, sequence = seq_obj.id, str(seq_obj.seq)

        # remove non-sequence character suffixes if they exist
        if not sequence[-1].isalpha() and sequence[-1] != '-':
            sequence = sequence[:-1]

        try:
            date_string = identifier.split('|')[2]
            date = np.datetime64(date_string)
            sequence_date_dict[sequence].append(date)
        except:
            pass

        sequence_count_dict[sequence] += 1
        sequence_length_dict[len(sequence)] += 1
        number_of_sequences += 1

    # sort sequences in decreasing order of frequency of occurrence
    sequence_count_dict = dict(sorted(sequence_count_dict.items(), key=lambda item: item[1], reverse=True))

    # compute median date of each sequence
    print('Computing median date for each sequence...')
    for sequence, date_list in tqdm(sequence_date_dict.items()):
        median_date = pd.Series(date_list, dtype='datetime64[ns]').quantile(0.5, interpolation="midpoint")
        sequence_date_dict[sequence] = median_date

    with open(os.path.join(data_directory, outfile), "w") as f:
        for seq, count in sequence_count_dict.items():
            seq_date = sequence_date_dict[seq]
            if seq_date is not None:
                print(f'>{count}|{seq_date}', file=f)
            else:
                print(f'>{count}', file=f)
            print(seq, file=f)

    print(f'Processed {infile}:')
    print(f'Found {number_of_sequences} sequences,')
    print(f'of which {len(sequence_count_dict)} are unique sequences.')

    return sequence_count_dict, sequence_length_dict


def create_variant_sequences_dict(sequences_src):
    """
    Creates a dictionary withs keys the variant names and values the reference spike protein sequences.
    """
    names = []
    for entry in os.listdir(sequences_src):  # Read all sequences
        if os.path.isfile(os.path.join(sequences_src, entry)):
            names.append(entry)

    consensus_dict = {}
    for name in names:
        if name[:-6] not in consensus_dict.keys():
            for fasta in SeqIO.parse(os.path.join(sequences_src, name), "fasta"):
                fasta_name, sequence = fasta.id, str(fasta.seq)
            consensus_dict[name[:-6]] = sequence

    return consensus_dict


def label_fasta_file_sequences_with_closest_variant(infile, outfile, path_to_consensus_sequences, data_directory):
    """
    Reads an unlabeled fasta file, finds the closest known variant to each sequence and
    generates an output file with each sequence labelled.
    """

    variant_sequences_dict = create_variant_sequences_dict(path_to_consensus_sequences)

    variant_similarity = {}
    unlabeled_fasta_file = open(os.path.join(data_directory, infile), "r")

    with open(os.path.join(data_directory, outfile), 'w') as labeled_fasta_file:
        while True:
            frequency = unlabeled_fasta_file.readline()
            sequence = unlabeled_fasta_file.readline()
            if not sequence: break  # EOF

            for variant, ref_seq in variant_sequences_dict.items():
                alignment_score = pairwise2.align.globalxx(ref_seq, sequence, score_only=True)
                variant_similarity[variant] = alignment_score

            line1 = frequency.strip() + "|" + max(variant_similarity, key=variant_similarity.get) + "\n"
            line2 = sequence
            labeled_fasta_file.writelines([line1, line2])


def reduce_and_align_sequences(infile: str,
                               outfile: str,
                               reduction_factor: int,
                               length_cutoff=1200,
                               invalid_amino_acids_cutoff=1,
                               data_directory=data_dir,
                               path_to_muscle_executable=path_to_muscle_executable):
    """
    Removes incomplete sequences from a fasta database, then downsamples, pools identical sequences, labels them with
    the closest matching covid spike protein variant and finally aligns them using MUSCLE.

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
    print('Labelling with closest variant...')
    label_fasta_file_sequences_with_closest_variant(infile=infile + '.cleaned.downsampled.unique',
                                                    outfile=infile + '.cleaned.downsampled.unique.labeled',
                                                    data_directory=data_dir,
                                                    path_to_consensus_sequences=path_to_consensus_sequences)

    # 5.
    print('To complete the alignment step run this command in the terminal:')
    muscle_command = MuscleCommandline(path_to_muscle_executable,
                                       input=data_directory + infile + '.cleaned.downsampled.unique.labeled',
                                       out=data_directory + outfile)
    print(muscle_command)


if __name__ == "__main__":
    # reduce_and_align_sequences(infile='spikeprot0112.fasta',
    #                            outfile='all_database_aligned.afa',
    #                            reduction_factor=1,
    #                            length_cutoff=1200,
    #                            invalid_amino_acids_cutoff=1)

    # combine_two_databases('./data/spike_protein_sequences/1_in_500_cleaned_aligned.afa',
    #                       'natural', './data/spike_protein_sequences/generated3.fasta', 'synthetic', 'combined.fasta')

    # print(data_dir)
    # muscle_command = MuscleCommandline(path_to_muscle_executable,
    #                   input='/home/dominic/PycharmProjects/VAE-SpikeProtein-Generation/data/spike_protein_sequences/spikeprot0112.fasta.cleaned.downsampled.unique.labeled',
    #                   out='/home/dominic/PycharmProjects/VAE-SpikeProtein-Generation/data/spike_protein_sequences/aligned.afa')
    # print(muscle_command)
    # remove_id_label_from_fasta_database(database='./data/spike_protein_sequences/spikeprot0112.fasta.cleaned.downsampled.unique.labeled',
    #                                     id_position=1,
    #                                     outfile='./data/spike_protein_sequences/spikeprot0112.fasta.cleaned.downsampled.unique.labeled.stripped')

    # partition_fasta_database_into_chunks('./data/spike_protein_sequences/spikeprot0112.fasta.cleaned.downsampled.unique.labeled.stripped',
    #                                      chunk_size=1500,
    #                                      outfile_prefix='./data/spike_protein_sequences/partitioned_database/spikeprot')

    # muscle_command = MuscleCommandline(path_to_muscle_executable,
    #                   input='path_to_infile',
    #                   out='path_to_outfile')
    # print(muscle_command)

    check_all_sequences_have_same_length('./data/spike_protein_sequences/1_in_50_cleaned.fasta')
