from Bio import SeqIO, pairwise2
from fasta_preprocessing_tools import combine_two_databases


def load_database(fasta_database):
    db = SeqIO.parse(fasta_database, 'fasta')

    return db


def get_conserved_regions_mask(database):

    for index, seq in enumerate(database):
        if index == 0:
            conserved_region_mask = [True]*len(seq.seq)
            reference_sequence = seq.seq

        else:
            for residue_index, residue in enumerate(seq.seq):
                if residue == reference_sequence[residue_index]:
                    pass
                else:
                    conserved_region_mask[residue_index] = False

    return reference_sequence, conserved_region_mask


def filter_by_conserved_regions(database, reference_sequence, conserved_region_mask):
    filtered_sequences = []
    for seq in database:
        for residue_index, residue in enumerate(seq.seq):
            if conserved_region_mask[residue_index] is True:
                if reference_sequence[residue_index] == residue:
                    pass
                else:
                    break
            if residue_index == len(seq.seq) - 1:
                filtered_sequences.append(seq)

    return filtered_sequences


if __name__ == "__main__":

    # Merge Generated Sequences into a Single File
    database1 = '../data/FcVAE_generated_high.fasta'
    database2 = '../data/FcVAE_generated_intermediate.fasta'
    database3 = '../data/FcVAE_generated_low.fasta'

    combine_two_databases(database1,
                          database2,
                          variant_database1='high',
                          variant_database2='intermediate',
                          outfile='../data/high_intermediate.fasta')

    combine_two_databases('../data/high_intermediate.fasta',
                          database3,
                          variant_database2='low',
                          outfile='../data/high_intermediate_low.fasta')

    # merge natural and synthetic databases
    combine_two_databases(database1='../data/high_intermediate_low.fasta',
                          database2='../data/spike_protein_sequences/1_in_500_cleaned_aligned.afa',
                          variant_database2='natural',
                          outfile='../data/merged_generated_and_natural.fasta')

    # muscle align the merged file


    # Load in generated sequences
    # generated_sequences_database = load_database('../data/high_intermediate_low.fasta')
    # natural_sequences_database = load_database('../data/spike_protein_sequences/1_in_500_cleaned_aligned.afa')
    #
    # # get list of conserved regions
    # ref_seq, conserved_region_mask = get_conserved_regions_mask(natural_sequences_database)
    # print(ref_seq)
    # print(conserved_region_mask)

    #filter by conserved regions
    #filtered_seqs = filter_by_conserved_regions(generated_sequences_database, ref_seq, conserved_region_mask)

    #print(len(filtered_seqs))