from Bio import pairwise2
import os
from Bio import SeqIO

cons_sequences_src = r"../data/spike_protein_sequences/consensus_sequences/"
sequences_src = r"C:\cdt_data\CovidProject\1_in_50_cleaned.fasta"

sequences_with_variant_src = r"C:\cdt_data\CovidProject\1_in_50_cleaned_variant.fasta"

def consensus_sequences_to_dict(sequences_src):
    names =[]
    for entry in os.listdir(sequences_src): #Read all sequences
        if os.path.isfile(os.path.join(sequences_src, entry)):
            names.append(entry)

    consensus_dict={}
    for name in names:
        if name[:-6] not in consensus_dict.keys():
            for fasta in SeqIO.parse(sequences_src + name, "fasta"):
                fasta_name, sequence = fasta.id, str(fasta.seq)
            consensus_dict[name[:-6]] = sequence
    return consensus_dict

consensus_dict = consensus_sequences_to_dict(cons_sequences_src)
print(consensus_dict["gamma"])

def sequences_to_id(file, consensus_dict):

        f = open(file, "r")
        while True:
            frequency = f.readline()
            sequence = f.readline()
            if not sequence: break  # EOF

            for variant, ref_seq  in consensus_dict.items():
                alignments = pairwise2.align.globalxx(ref_seq, sequence)

        with open(sequences_with_variant_src, 'a') as outfile:
            line1 = "PS5 Restock India \n"
            line2 = "Xbox Series X Restock India \n"
            line3 = "Nintendo Switch Restock India"

            outfile.writelines([line1, line2, line3])

sequences_to_id(sequences_src)



#print(names)
def calculate_similarity_scores():
    alignments = pairwise2.align.globalxx("ACCGT", "ACG")