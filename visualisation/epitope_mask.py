import numpy as np
from Bio import SeqIO
from DSSPparser import parseDSSP
import pandas as pd
from collections import defaultdict
from Bio import pairwise2

#TODO: do we even need this function
def convert_fasta_file_to_df(fasta_file):
    """
    Reads a fasta file containing a single sequence and coverts that sequence to a pandas data frame.
    """
    number_of_sequences_found = 0
    sequence_length = None
    sequence = None
    for sequence_entry in SeqIO.parse(fasta_file, "fasta"):
        name, sequence = sequence_entry.id, str(sequence_entry.seq)
        number_of_sequences_found += 1
        sequence_length = len(sequence)

    if number_of_sequences_found != 1:
        raise ValueError('fasta file should have only a single sequence')

    # Create an ID list for the sequences
    pos = [str(i) for i in list(range(1, sequence_length+1))]

    seq_df = pd.DataFrame(list(zip(sequence, pos)), columns=['amino_acid', 'position'])

    return seq_df


def extract_solvent_accessibility_from_dssp_file(dssp_file):
    '''
    TODO: description
    '''
    dssp_parser = parseDSSP(dssp_file)
    dssp_parser.parse()
    dssp_parsed_dict = dssp_parser.dictTodataframe()
    dssp_parsed_dict = dssp_parsed_dict[['resnum', 'inscode', 'chain', 'aa', 'acc']]

    solvent_accessibility_of_residue = defaultdict(lambda: [])
    residue_number_to_amino_acid = defaultdict(lambda: '')
    for entry in dssp_parsed_dict.iterrows():
        # TODO: comment why [1]
        amino_acid = entry[1].aa
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


def align_dssp_to_fasta_and_return_epitope_mask(fasta_df, dssp_df):
    """
    TODO: description
    """
    fasta_sequence = pandas_series_to_string(fasta_df.amino_acid)
    dssp_sequence = pandas_series_to_string(dssp_df.amino_acid)

    alignments = pairwise2.align.globalxx(fasta_sequence, dssp_sequence)

    #TODO: check
    best_alignment = alignments[0]

    aligned_fasta_sequence = best_alignment.seqA
    aligned_dssp_sequence = best_alignment.seqB

    #TODO: rename
    pointer_in_dssp_sequence = 0
    fasta_solvent_accessibility_vector = []

    for index, dssp_amino_acid in enumerate(aligned_dssp_sequence):
        if aligned_fasta_sequence[index] == "-":
            continue
        else:
            if dssp_amino_acid != "-":
                solvent_accessibility_of_amino_acid = dssp_df.solvent_accessibility[pointer_in_dssp_sequence]
                fasta_solvent_accessibility_vector.append(solvent_accessibility_of_amino_acid)
                pointer_in_dssp_sequence += 1
            else:
                fasta_solvent_accessibility_vector.append(None)

    return fasta_solvent_accessibility_vector


def extract_solvent_accessibility_vector_from_fasta_and_dssp_file(fasta_file, dssp_file, data_directory):
    """
    TODO: description
    """
    fasta_df = convert_fasta_file_to_df(data_directory+fasta_file)
    dssp_df = extract_solvent_accessibility_from_dssp_file(data_directory + dssp_file)

    return align_dssp_to_fasta_and_return_epitope_mask(fasta_df, dssp_df)


if __name__ == "__main__":
    data_dir = "../data/spike_protein_pdb/"
    spike_structure_src = r"../data/spike_protein_pdb/7n1u.pdb"
    spike_sequence_src = r"../data/spike_protein_pdb/rcsb_pdb_7N1U.fasta"
    spike_dssp_src = r"../data/spike_protein_pdb/7n1u.dssp"

    SAV = extract_solvent_accessibility_vector_from_fasta_and_dssp_file(fasta_file='rcsb_pdb_7N1U.fasta',
                                                                        dssp_file='7n1u.dssp',
                                                                        data_directory=data_dir)

    print(SAV)

