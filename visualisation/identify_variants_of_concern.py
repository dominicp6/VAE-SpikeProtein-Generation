from Bio import pairwise2
import os
from Bio import SeqIO


def consensus_sequences_to_dict(sequences_src):
    names = []
    for entry in os.listdir(sequences_src):  # Read all sequences
        if os.path.isfile(os.path.join(sequences_src, entry)):
            names.append(entry)

    consensus_dict = {}
    for name in names:
        if name[:-6] not in consensus_dict.keys():
            for fasta in SeqIO.parse(sequences_src + name, "fasta"):
                fasta_name, sequence = fasta.id, str(fasta.seq)
            consensus_dict[name[:-6]] = sequence
    return consensus_dict


def sequences_to_id(file, consensus_dict):
    variant_similarity = {}
    f = open(file, "r")
    while True:
        frequency = f.readline()
        # print(frequency)
        sequence = f.readline()
        if not sequence: break  # EOF

        for variant, ref_seq in consensus_dict.items():
            alignment_score = pairwise2.align.globalxx(ref_seq, sequence, score_only=True)
            variant_similarity[variant] = alignment_score

        with open(sequences_with_variant_src, 'a') as outfile:
            line1 = frequency.strip() + "," + max(variant_similarity, key=variant_similarity.get) + "\n"
            line2 = sequence
            outfile.writelines([line1, line2])


if __name__ == "__main__":
    cons_sequences_src = r"../data/spike_protein_sequences/consensus_sequences/"
    sequences_src = r"C:\cdt_data\CovidProject\1_in_50_cleaned.fasta"
    # sequences_with_variant_src = r"C:\cdt_data\CovidProject\1_in_500_cleaned_variant.fasta"
    sequences_with_variant_src = sequences_src.split('.')[0] + "_variant.fasta"

    consensus_sequences_dict = consensus_sequences_to_dict(cons_sequences_src)

    sequences_to_id(sequences_src, consensus_sequences_dict)
