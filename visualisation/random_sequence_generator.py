"""
Generates "random" spike protein sequences using different methods.

NB: 1_in_500_cleaned_unique.afa has 8880 sequences, of which 2274 are unique
"""

import numpy as np
from Bio import SeqIO
from nltk.lm.preprocessing import padded_everygram_pipeline
from nltk.lm import MLE
from tqdm import tqdm


sequence_length = 1000

valid_residue_types = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                       'R', 'S', 'T', 'V', 'W', 'Y', '-']


def generate_random_ngram_sequences(fasta_file, N, L, outfile, n=3):
    lm = MLE(n)
    fasta = SeqIO.parse(fasta_file, 'fasta')

    print('Reading data file')
    sequence_data = []
    for sequence in fasta:
        sequence_data.append([letter for letter in sequence.seq])

    print('Training model')
    train, vocab = padded_everygram_pipeline(n, sequence_data)
    lm.fit(train, vocab)

    print('Generating sequences')
    with open(outfile, "w") as f:
        for _ in tqdm(range(N)):
            while True:
                seq = lm.generate(L-(n-1), text_seed='MFVFLVLLPLVSSQCVN'[0:n-1])
                seq.insert(0, 'F')
                seq.insert(0, 'M')
                seq = "".join(seq)
                if "</s>" not in seq:      # only generate sequences that dont have a stop symbol in them
                    print(">", file=f)
                    print(seq, file=f)
                    break                  # break the while loop

    return lm


def generate_completely_random_sequences(N, L, outfile):
    with open(outfile, "w") as f:
        for _ in range(N):
            seq = "".join(np.random.choice(valid_residue_types, size=L))
            print(">", file=f)
            print(seq, file=f)


if __name__ == "__main__":
    #generate_sequences(100, 1282, outfile="completely_random_sequences")
    lm = generate_random_ngram_sequences("1_in_500_cleaned_aligned.afa", 8880, 1282, n=18, outfile='random_18gram_sequences')
    from fasta_preprocessing_tools import reduce_to_unique_sequences
    reduce_to_unique_sequences(infile="random_18gram_sequences",
                               outfile="random_18gram_sequences.unique",
                               data_directory='.')

