from Bio import SeqIO
from DSSPparser import parseDSSP
import pandas as pd
from collections import defaultdict
from Bio import pairwise2


def convert_fasta_file_to_df(fasta_src):
    number_of_sequences_found = 0
    sequence_length = None
    sequence = None
    for sequence_entry in SeqIO.parse(fasta_src, "fasta"):
        name, sequence = sequence_entry.id, str(sequence_entry.seq)
        number_of_sequences_found += 1
        sequence_length = len(sequence)

    if number_of_sequences_found != 1:
        raise ValueError('fasta file should have only a single sequence')

    # Create an ID list for the sequences
    pos = [str(i) for i in list(range(1, sequence_length+1))]

    seq_df = pd.DataFrame(list(zip(sequence, pos)), columns=['residue', 'position'])

    return seq_df


def extract_SASA_from_dssp(dssp_src):
    '''
    #print(pddict.columns)
    #Index(['resnum', 'inscode', 'chain', 'aa', 'struct', 'structdetails', 'bp1',
    #       'bp2', 'acc', 'h_nho1', 'h_ohn1', 'h_nho2', 'h_ohn2', 'tco', 'kappa',
    #   'alpha', 'phi', 'psi', 'xca', 'yca', 'zca', 'rcsb_given_chain',
    #        'author_given_chain'],      dtype='object'
    :param dssp_src:
    :return:
    '''
    dssp_parser = parseDSSP(dssp_src)
    dssp_parser.parse()
    dssp_parsed_dict = dssp_parser.dictTodataframe()
    sasa_from_dssp = dssp_parsed_dict[['resnum', 'inscode', 'chain', 'aa', 'acc']]
    #sasa_from_dssp = dssp_parsed_dict[['inscode', 'acc']]
    #sasa_from_dssp = sasa_from_dssp.groupby(by = ['inscode']).median()

    solvent_accessibility_dictionary = defaultdict(lambda: [])
    for entry in sasa_from_dssp.iterrows():
        amino_acid_residue_number = entry[1].inscode
        amino_acid = entry[1].aa
        solvent_accessibility = entry[1].acc
        solvent_accessibility_dictionary[amino_acid_residue_number].append(solvent_accessibility)

    averaged_solvent_accessibility = dict()
    for residue_number, solvent_accessibilities in solvent_accessibility_dictionary.items():
        averaged_solvent_accessibility[residue_number] = np.median(solvent_accessibilities)

    return averaged_solvent_accessibility


def align_dssp_to_fasta_and_return_epitope_mask(fasta_df, dssp_df):
    fasta_sequence = fasta_df.residue
    dssp_sequence = dssp_df.aa

    print(fasta_sequence)
    print(dssp_sequence)



    alignments = pairwise2.align.globalxx(fasta_sequence, dssp_sequence)

    best_alignment = alignments[0]

    print(best_alignment)
    print(best_alignment.__vars__)




if __name__ == "__main__":
    spike_structure_src = r"../data/spike_protein_pdb/7n1u.pdb"
    spike_sequence_src = r"../data/spike_protein_pdb/rcsb_pdb_7N1U.fasta"
    spike_dssp_src = r"../data/spike_protein_pdb/7n1u.dssp"

    # 1. Read and covert pdb to dssp
    # Terminal command:  mkdssp -i original_covid.pdb -o original_covid.dssp

    # 2. Extract solvent accessibility from dssp
    print('DSSP')
    dssp_df = extract_SASA_from_dssp(spike_dssp_src)
    dssp_df.to_csv(r"../data/spike_protein_pdb/dssp_test.csv", index=False, na_rep='NULL')
    print(dssp_df.dtypes)

    # 3. Align fasta with dssp solvent accessibility, output epitope mask
    print('FASTA')
    fasta_df = convert_fasta_file_to_df(spike_sequence_src)
    fasta_df.to_csv(r"../data/spike_protein_pdb/fasta_test.csv", index=False, na_rep='NULL')
    print(fasta_df)

    align_dssp_to_fasta_and_return_epitope_mask(fasta_df, dssp_df)
    # df = fasta_df.merge(dssp_df, how="left", left_on='position', right_on='inscode')
    # print(df.head(20))
    #
    # epitope_mask = df[['sequence', 'acc']]
    #
    # df.to_csv(r"../data/spike_protein_pdb/epitope_mask.csv", index=False, na_rep='NULL')


