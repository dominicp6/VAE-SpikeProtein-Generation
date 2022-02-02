spike_structure_src = r"../data/spike_protein_pdb/7n1u.pdb"
spike_sequence_src = r"../data/spike_protein_pdb/rcsb_pdb_7N1U.fasta"
spike_dssp_src = r"../data/spike_protein_pdb/7n1u.dssp"

from Bio.PDB import *
from Bio import SeqIO
from DSSPparser import parseDSSP
from biopandas.pdb import PandasPdb
import pandas as pd

# from Bio import pairwise2
#     alignments = pairwise2.align.globalxx("ACCGT", "ACG")



def fasta_to_df(fasta_src):
    for fasta in SeqIO.parse(fasta_src, "fasta"):
        name, sequence = fasta.id, str(fasta.seq)
    pos = [str(i) for i in list(range(1, len(sequence)+1))]
    #print(pos)
    seq_df = pd.DataFrame(list(zip(sequence, pos)), columns= ['sequence','position'])
    return seq_df

def dssp_to_SASA(dssp_src):
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
    pddict = dssp_parser.dictTodataframe()
    sasa_from_dssp = pddict[['resnum', 'inscode', 'chain', 'aa', 'acc']]
    #sasa_from_dssp = pddict[['inscode', 'acc']]
    #sasa_from_dssp = sasa_from_dssp.groupby(by = ['inscode'])[['acc']].median()
    return sasa_from_dssp

dssp_df = dssp_to_SASA(spike_dssp_src)

fasta_df = fasta_to_df(spike_sequence_src)

dssp_df.to_csv(r"../data/spike_protein_pdb/dssp_test.csv", index=False, na_rep='NULL')
fasta_df.to_csv(r"../data/spike_protein_pdb/fasta_test.csv", index=False, na_rep='NULL')
print(fasta_df)

print(dssp_df.dtypes)

df = fasta_df.merge(dssp_df, how="left", left_on='position', right_on='inscode')
print(df.head(20))

epitope_mask = df[['sequence', 'acc']]

df.to_csv(r"../data/spike_protein_pdb/epitope_mask.csv", index=False, na_rep='NULL')

print(epitope_mask.head(20))


#for model in structure:
#    for chain in model:
        #print(chain)
 #       for residue in chain:
          # print(residue)
            #for atom in residue:
                #print(atom)


#Terminal command:  mkdssp -i original_covid.pdb -o original_covid.dssp

# Read in dssp file and construct solvent accessibility vector

def extract_solvent_accessibility_vector_from_dssp(sequence_src, structure_src):

    structure = parser.get_structure('X', structure_src)
    model = structure[0]
    #dssp = DSSP(model, spike_structure_src)
    #print(vars(dssp))

    resolution = structure.header['resolution']
    keywords = structure.header['keywords']
    print(resolution, keywords)


extract_solvent_accessibility_vector_from_dssp(spike_sequence_src, spike_structure_src)




def pdb_to_residues(pdb_src):
    ppdb = PandasPdb()
    pdb_df = ppdb.read_pdb(pdb_src)
    print(pdb_df.df['ATOM'].head(3))

    print(pdb_df)
    residues_from_pdb = pdb_df[[""]]

    return residues_from_pdb



print(pdb_to_residues(spike_structure_src))

parser = PDBParser()

structure = parser.get_structure('X',spike_structure_src)
#io = PDBIO()
print(structure['REMARK'])

