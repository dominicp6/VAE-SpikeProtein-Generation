import numpy as np
from Bio import SeqIO
from DSSPparser import parseDSSP
import pandas as pd
from collections import defaultdict
from Bio import pairwise2


def get_fasta_sequence_from_file(fasta_file):
    """
    Reads a fasta file containing a single sequence and returns the sequence string.
    """
    number_of_sequences_found = 0
    sequence = None
    for sequence_entry in SeqIO.parse(fasta_file, "fasta"):
        _, sequence = sequence_entry.id, str(sequence_entry.seq)
        number_of_sequences_found += 1

    if number_of_sequences_found != 1:
        raise ValueError('fasta file should have only a single sequence')

    return sequence


def extract_solvent_accessibility_from_dssp_file(dssp_file):
    """
    Parses a dssp file to find the solvent accessibility of each residue.
    """
    dssp_parser = parseDSSP(dssp_file)
    dssp_parser.parse()
    dssp_parsed_dict = dssp_parser.dictTodataframe()
    dssp_parsed_dict = dssp_parsed_dict[['resnum', 'inscode', 'chain', 'aa', 'acc']]

    solvent_accessibility_of_residue = defaultdict(lambda: [])
    residue_number_to_amino_acid = defaultdict(lambda: '')
    for entry in dssp_parsed_dict.iterrows():
        amino_acid = entry[1].aa  # position [0] is ID, position [1] is row data

        # the ! character is used to mark sequence breaks in dssp files; we ignore it
        if amino_acid == "!":
            continue

        # In dssp files CYS amino acids are labelled with lower case letters
        if amino_acid.islower():
            # relabel lower case letter with one-letter code for CYS
            amino_acid = 'C'

        amino_acid_residue_number = int(entry[1].inscode)

        residue_number_to_amino_acid[amino_acid_residue_number] = amino_acid
        solvent_accessibility = entry[1].acc
        solvent_accessibility_of_residue[amino_acid_residue_number].append(int(solvent_accessibility))

    averaged_solvent_accessibility_of_residue = dict()
    for residue_number, solvent_accessibilities in solvent_accessibility_of_residue.items():
        averaged_solvent_accessibility_of_residue[residue_number] = np.median(solvent_accessibilities)

    dssp_dictionary = {'residue_number': list(averaged_solvent_accessibility_of_residue.keys()),
                       'solvent_accessibility': list(averaged_solvent_accessibility_of_residue.values()),
                       'amino_acid': [residue_number_to_amino_acid[residue]
                                      for residue in averaged_solvent_accessibility_of_residue.keys()]}

    return pd.DataFrame(dssp_dictionary)


def pandas_series_to_string(series):
    """
    Concatenates the entries of a pandas series into a string.
    """
    string = ""
    for entry in series:
        string += str(entry)

    return string


def align_sequences(sequence1, sequence2):
    alignments = pairwise2.align.globalxx(sequence1, sequence2)

    best_alignment = alignments[0]

    sequence1_aligned = best_alignment.seqA
    sequence2_aligned = best_alignment.seqB

    return sequence1_aligned, sequence2_aligned


def get_solvent_accessibility_vector_from_fasta_and_dssp(fasta_file, dssp_file, data_directory):
    """
    Pairwise sequence aligns a reference fasta sequence to a dssp sequence and returns
    a solvent accessibility vector for the fasta sequence.
    """
    fasta_sequence = get_fasta_sequence_from_file(data_directory+fasta_file)
    dssp_solvent_accessibility = extract_solvent_accessibility_from_dssp_file(data_directory+dssp_file)
    dssp_sequence = pandas_series_to_string(dssp_solvent_accessibility.amino_acid)

    aligned_fasta_sequence, aligned_dssp_sequence = align_sequences(fasta_sequence, dssp_sequence)

    position_in_dssp_sequence = 0
    fasta_solvent_accessibility_vector = []

    for index, dssp_amino_acid in enumerate(aligned_dssp_sequence):
        if aligned_fasta_sequence[index] == "-":
            continue
        else:
            if dssp_amino_acid != "-":
                solvent_accessibility_of_amino_acid = dssp_solvent_accessibility.solvent_accessibility[position_in_dssp_sequence]
                fasta_solvent_accessibility_vector.append(solvent_accessibility_of_amino_acid)
                position_in_dssp_sequence += 1
            else:
                fasta_solvent_accessibility_vector.append(0)

    return fasta_solvent_accessibility_vector


if __name__ == "__main__":
    data_dir = "../data/spike_protein_pdb/"
    spike_structure_src = r"../data/spike_protein_pdb/7n1u.pdb"
    spike_sequence_src = r"../data/spike_protein_pdb/rcsb_pdb_7N1U.fasta"
    spike_dssp_src = r"../data/spike_protein_pdb/7n1u.dssp"

    SAV = get_solvent_accessibility_vector_from_fasta_and_dssp(fasta_file='rcsb_pdb_7N1U.fasta',
                                                               dssp_file='7n1u.dssp',
                                                               data_directory=data_dir)

    print(SAV)
