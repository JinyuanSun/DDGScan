#!/usr/bin/env python
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)

def main(structure_id, filename, chain, resnum):

    structure = parser.get_structure(structure_id, filename)

    model = structure[0]
    chain = model[chain]
    residue = chain[resnum]




if __name__ == '__main__':
    structure_id = "1LVE"
    filename = "../1LVE_modfixed.pdb"
    chain = 'A'
