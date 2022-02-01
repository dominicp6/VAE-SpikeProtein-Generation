spike_structure_src = r"../data/spike_protein_pdb/original_covid.pdb"

from Bio.PDB import *

parser = PDBParser()

io = PDBIO()

structure = parser.get_structure('X', spike_structure_src)
model = structure[0]
dssp = DSSP(model, spike_structure_src)
print(vars(dssp))

#Terminal command:  mkdssp -i original_covid.pdb -o original_covid.dssp


# Read in dssp file and construct solvent accessibility vector

def extract_solvent_accessibility_vector_from_dssp(infile):
