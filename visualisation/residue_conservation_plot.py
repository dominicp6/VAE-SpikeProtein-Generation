"""
A script to visualise how much each residue of the Spike protein is conserved across the population.

To compute spike protein variation at each residue one method is to follow Sikora et al
(except I think that we can improve on it! They fixed their spike proteins to be of a given length!):

Define a conservation score to be 1 minus the entropy of the observed amino acid distribution to the maximum possible
entropy at a given position: Score(i) = 1-\left(\frac{-\sum_k p_k(i)\log p_k(i)}{\log(21)}\right)

Note log(21) not log(20) because we include the empty code "-" as a possibility since we are dealing with aligned
sequences.
"""

from Bio import SeqIO
from collections import defaultdict
import os
import numpy as np
import matplotlib.pyplot as plt

script_dir = os.path.dirname(os.path.realpath(__file__))  # path to this file (DO NOT CHANGE)

# relative path of datasets
data_dir = script_dir + '/data/spike_protein_sequences/'

fasta_file = "../data/spike_protein_sequences/1_in_500_cleaned_aligned.afa"
fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')

residue_distribution = defaultdict(lambda: {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "Q": 0, "E": 0, "G": 0, "H": 0, "I": 0, "L": 0,
                         "K": 0, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0, "-": 0})
for fasta in fasta_sequences:
    id_label, sequence = fasta.id, str(fasta.seq)
    sequence_count = int(id_label.split('|')[0])
    for position, letter in enumerate(sequence):
        residue_distribution[position][letter] += sequence_count


def calculate_entropy_of_count_distribution(count_distribution):
    normalised_distribution = count_distribution/sum(count_distribution)
    distribution_entropy = sum([probability * np.log(probability)
                                for probability in normalised_distribution if probability > 0])

    return distribution_entropy


def calculate_conservation_score(entropy, maximum_entropy):
    return 1 - entropy/maximum_entropy


conservation_vector = np.zeros(len(residue_distribution))
for position, count_dictionary in residue_distribution.items():
    entropy_of_position = calculate_entropy_of_count_distribution(np.fromiter(count_dictionary.values(), dtype=np.int64))
    position_conservation_score = calculate_conservation_score(entropy_of_position, maximum_entropy=-np.log(21))
    conservation_vector[position] = position_conservation_score


print(len(conservation_vector))
print(np.count_nonzero(conservation_vector == 1))
plt.plot(1-(1-conservation_vector)/np.max(1-conservation_vector), 'k', linewidth=0.75)
plt.xlabel('Residue Number', fontsize=20)
plt.ylabel('Conservation Score', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.margins(x=0, y=0)
plt.ylim((0.00, 1.00))
fig = plt.gcf()
fig.set_size_inches(16.5, 10.5)
fig.set_dpi(100)
plt.show()




