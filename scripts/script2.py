'''
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
p = PDBParser()
structure = p.get_structure("APLHA", "./ranked_0.pdb")
model = structure[0]
dssp = DSSP(model, "./ranked_0.pdb")

a_key = list(dssp.keys())[2]

dssp[a_key]


'''
from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename='ranked_0.pdb')
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)

PDBFile.writeFile(fixer.topology, fixer.positions, open('output.pdb', 'w'))