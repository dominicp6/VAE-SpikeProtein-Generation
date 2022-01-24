# Adapted from https://github.com/xyjing-works/SequenceEncoding/blob/master/SequenceEncoding.py

import json


class SequenceEncoding:
    encoding_types = ['One_hot', 'One_hot_6_bit', 'Binary_5_bit', 'Hydrophobicity_matrix',
                      'Meiler_parameters', 'Acthely_factors', 'PAM250', 'BLOSUM62', 'Miyazawa_energies',
                      'Micheletti_potentials', 'AESNN3', 'ANN4D', 'ProtVec']
    residue_types = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
                     'X']

    def __init__(self, encoding_type="One_hot"):
        if encoding_type not in SequenceEncoding.encoding_types:
            raise Exception("Encoding type \'%s\' not found" % encoding_type)
        self.encoding_type = encoding_type

    def get_ProtVec_encoding(self, ProtVec, seq, overlap=True):
        if overlap:
            encodings = []
            for i in range(len(seq) - 2):
                encodings.append({seq[i:i + 3]: ProtVec[seq[i:i + 3]]}) if ProtVec.__contains__(
                    seq[i:i + 3]) else encodings.append({seq[i:i + 3]: ProtVec["<unk>"]})
        else:
            encodings_1, encodings_2, encodings_3 = [], [], []
            for i in range(0, len(seq), 3):
                if i + 3 <= len(seq):
                    encodings_1.append({seq[i:i + 3]: ProtVec[seq[i:i + 3]]}) if ProtVec.__contains__(
                        seq[i:i + 3]) else encodings_1.append({seq[i:i + 3]: ProtVec["<unk>"]})
                if i + 4 <= len(seq):
                    encodings_2.append({seq[i + 1:i + 4]: ProtVec[seq[i + 1:i + 4]]}) if ProtVec.__contains__(
                        seq[i + 1:i + 4]) else encodings_2.append({seq[i + 1:i + 4]: ProtVec["<unk>"]})
                if i + 5 <= len(seq):
                    encodings_3.append({seq[i + 2:i + 5]: ProtVec[seq[i + 2:i + 5]]}) if ProtVec.__contains__(
                        seq[i + 2:i + 5]) else encodings_3.append({seq[i + 2:i + 5]: ProtVec["<unk>"]})

            encodings = [encodings_1, encodings_2, encodings_3]
        return encodings

    def get_encoding(self, seq, overlap=True):
        seq = seq.upper()
        with open("data/encodings/%s.json" % self.encoding_type, 'r') as load_f:
            encoding = json.load(load_f)
        encoding_data = []
        if self.encoding_type == "ProtVec":
            encoding_data = self.get_ProtVec_encoding(encoding, seq, overlap)
        else:
            for res in seq:
                if res not in SequenceEncoding.residue_types:
                    res = "X"
                encoding_data.append({res: encoding[res]})

        return encoding_data

    def encode(self, seq):
        return [digit for dictionary in self.get_encoding(seq) for amino_acid, encoding_vector in dictionary.items() for
                digit in encoding_vector]


if __name__ == "__main__":
    test_spike_protein_sequence = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAISGTNGTKRFDNPVLPFNDGVYFASTEKSNI" \
                                  "IRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDG" \
                                  "YFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLS" \
                                  "ETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNV" \
                                  "YADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQ" \
                                  "SYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIDDTTDAVRDPQTLEILDITPCS" \
                                  "FGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSHRRARSVAS" \
                                  "QSIIAYTMSLGAENSVAYSNNSIAIPINFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQ" \
                                  "IYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSG" \
                                  "WTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILAR" \
                                  "LDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAI" \
                                  "CHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTHNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNI" \
                                  "QKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"

    # example One-Hot encoding
    one_hot_encoder = SequenceEncoding(encoding_type='One_hot')
    encoded_sequence = one_hot_encoder.encode(test_spike_protein_sequence)
    print(encoded_sequence)

    # example PAM250 encoding
    PAM250_encoder = SequenceEncoding(encoding_type='PAM250')
    encoded_sequence = PAM250_encoder.encode(test_spike_protein_sequence)
    print(encoded_sequence)
