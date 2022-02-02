
from Bio import pairwise2
import os
from Bio import SeqIO




sequences_src = r"../data/spike_protein_sequences/consensus_sequences/"

names =[]
for entry in os.listdir(sequences_src): #Read all sequences
    if os.path.isfile(os.path.join(sequences_src, entry)):
        names.append(entry)

data_dict={}
for name in names:
    if name[:-6] not in data_dict.keys():
        for fasta in SeqIO.parse(fasta_src, "fasta"):
            fasta_name, sequence = fasta.id, str(fasta.seq)
            data_dict[name]=sequence
            print(data_dict)



print(names)
def calculate_similarity_scores():
    alignments = pairwise2.align.globalxx("ACCGT", "ACG")