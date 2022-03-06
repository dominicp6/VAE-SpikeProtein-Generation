import numpy as np
from tqdm import tqdm
from Bio import pairwise2, SeqIO



def find_closest_natural_sequence(synthetic_sequence, natural_sequences_database):
    minimum_alignment_score = np.float('inf')
    aligned_natural_seq = ''
    aligned_synthetic_seq = ''
    for natural_seq in tqdm(natural_sequences_database):
        alignment = pairwise2.align.globalxx(synthetic_sequence, natural_seq.seq)
        best_alignment = alignment[0]
        alignment_score = best_alignment.score

        if alignment_score < minimum_alignment_score:
            minimum_alignment_score = alignment_score
            aligned_synthetic_seq = best_alignment.seqA
            aligned_natural_seq = best_alignment.seqB
        else:
            pass

    return aligned_synthetic_seq, aligned_natural_seq, minimum_alignment_score



def generate_point_mutation_list(original_sequence, mutated_sequence):
    point_mutation_list = ""
    for res_index, (orig_res, mut_res) in enumerate(zip(original_sequence, mutated_sequence)):
        if mut_res == "-" or orig_res == "-":
            # potentially an insertion or deletion -ignore
            continue
        elif mut_res != orig_res:
            # a mutation - add it to the list
            point_mutation_list+=f"{orig_res}{res_index+1}{mut_res} \n"
        else:
            # no mutation - ignore
            continue

    return point_mutation_list



if __name__ == "__main__":
    synthetic_db = SeqIO.parse('../data/high_intermediate_low.fasta', 'fasta')
    natural_db = SeqIO.parse('../data/spike_protein_sequences/1_in_500_cleaned_aligned.afa', 'fasta')

    for seq in synthetic_db:
        aligned_synthetic, aligned_natural, minimum_alignment_score = find_closest_natural_sequence(seq.seq, natural_db)

        print(f'Found a closely aligned sequence with alignment score {minimum_alignment_score}.')
        print(f'nat: {aligned_natural}')
        print(f'syn: {aligned_synthetic}')
        print(generate_point_mutation_list(original_sequence=aligned_natural, mutated_sequence=aligned_synthetic))
        input()


