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
    pdbid = pdb[:4]
    os.system("mv %s %s_raw.pdb" % (pdbid, pdbid))
    # Retrieve structure from PDB.
    print('Retrieving %s from PDB...' % pdbid)
    fixer = PDBFixer(pdb)

    # Build a list of chains to remove.
    print('Removing all chains but %s' % chain_ids_to_keep)
    # all_chains = list(fixer.topology.chains())

    chain_id_list = [c.id for c in fixer.topology.chains()]
    chain_ids_to_remove = set(chain_id_list) - set(chain_ids_to_keep)
    fixer.removeChains(chainIds=chain_ids_to_remove)
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
    missing_atoms = fixer.missingAtoms
    # print(missing_atoms)

    need_atoms = {}
    for residue, atoms in missing_atoms.items():
        for atom in atoms:
            if atom.name in ["C", "N", "CA", "O"]:
                # keep mainchain residues only
                need_atoms[residue] = [atom]

    fixer.missingAtoms = need_atoms
    # print(fixer.missingAtoms)
    fixer.addMissingAtoms()

    # Remove heterogens.
    print('Removing heterogens...')
    fixer.removeHeterogens(keepWater=False)

    # Write PDB file.
    output_filename = '%s.pdb' % pdbid
    print('Writing PDB file to "%s"...' % output_filename)
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_filename, 'w'))
    return output_filename


if __name__ == '__main__':
    pdb = '4lvr_raw.pdb'  # PDB ID to retrieve
    chain_ids_to_keep = ['A']  # chains to keep
    autofix(pdb, chain_ids_to_keep)
