from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

fixer = PDBFixer('6JTT_A_prot.pdb')
#fixer.removeChains([1,2,3,4,5,7,9])
fixer.findMissingResidues()

# only add missing residues in the middle of the chain, do not add terminal ones
chains = list(fixer.topology.chains())
keys = fixer.missingResidues.keys()
missingResidues = dict()
for key in keys:
    chain = chains[key[0]]
    if not (key[1] == 0 or key[1] == len(list(chain.residues()))):
        missingResidues[key] = fixer.missingResidues[key]
fixer.missingResidues = missingResidues

fixer.findMissingAtoms()
fixer.addMissingAtoms()

PDBFile.writeFile(fixer.topology, fixer.positions, open('6JTT_fixed.pdb', 'w'))
