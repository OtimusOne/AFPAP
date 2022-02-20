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