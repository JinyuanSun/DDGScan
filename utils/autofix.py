#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/3/28 6:49 PM
# @File    : autofix.py
# @annotation    : Fix missing atoms in mainchain
# !/usr/bin/env python

"""

"""

import os

from openmm import app
from pdbfixer import PDBFixer


################################################################################
# ref:https://python.hotexamples.com/site/file?hash=0xc5901342b2c339b10661d6508\
# 4166c386ffc205460660dd91e9b1ce00555c106&fullName=openmmtools/data/prepare_pdb\
# .py&project=choderalab/openmmtools
################################################################################


def write_file(filename, contents):
    outfile = open(filename, 'w')
    outfile.write(contents)
    outfile.close()


################################################################################
# SET UP SYSTEM
################################################################################

def autofix(pdb, chain_ids_to_keep, risky=False):

    os.system("cp %s %s.raw" % (pdb, pdb))
    fixer = PDBFixer(filename=pdb)
    fixer_1 = PDBFixer(filename=pdb)
    fixer_2 = PDBFixer(filename=pdb)

    # Build a list of chains to remove.
    print('Removing all chains but %s' % chain_ids_to_keep)
    # all_chains = list(fixer.topology.chains())

    chain_id_list = [c.id for c in fixer.topology.chains()]
    chain_ids_to_remove = set(chain_id_list) - set(chain_ids_to_keep)
    fixer.removeChains(chainIds=chain_ids_to_remove)
    fixer_1.removeChains(chainIds=chain_ids_to_remove)
    fixer_2.removeChains(chainIds=chain_ids_to_remove)
    if risky:
        # Replace nonstandard residues.
        print('Replacing nonstandard residues...')
        print("It is risky to replace ncAA using pdbfixer!")
        # logging.WARNING("It is risky to replace ncAA using pdbfixer!")
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

    # Add missing atoms.
    print('Adding backbone missing atoms only!')

    fixer.findMissingResidues()

    # modeller = app.Modeller(self.topology, self.positions)
    # modeller.delete(toDelete)
    fixer.findMissingAtoms()
    print(f"MISSING ATOMS 000: {fixer.missingAtoms}")
    missing_atoms = fixer.missingAtoms
    # print(missing_atoms)

    need_atoms = {}
    for residue, atoms in missing_atoms.items():
        for atom in atoms:
            if atom.name in ["C", "N", "CA", "O"]:
                # keep mainchain residues only
                need_atoms[residue] = [atom]

    fixer_1.missingAtoms = need_atoms
    print(f"MISSING ATOMS 111: {fixer_1.missingAtoms}")
    print(f"MISSING terminals: =======\n {fixer.missingTerminals} ======\n")
    fixer_1.missingTerminals = fixer.missingTerminals
    fixer_1.missingResidues = need_atoms.keys()
    fixer_1.addMissingAtoms()

    # Remove heterogens.
    print('Removing heterogens...')
    fixer_1.removeHeterogens(keepWater=False)

    # Write PDB file.
    output_filename = pdb.replace(".pdb", "_fixed.pdb")
    print('Writing PDB file to "%s"...' % output_filename)
    app.PDBFile.writeFile(fixer_2.topology, fixer_2.positions, open(output_filename, 'w'))
    return output_filename


if __name__ == '__main__':
    pdb = '../test_bak/1NWW.pdb'  # PDB ID to retrieve
    chain_ids_to_keep = ['A']  # chains to keep
    autofix(pdb, chain_ids_to_keep)
