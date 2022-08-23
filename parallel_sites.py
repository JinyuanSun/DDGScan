#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/8/8 10:02 AM
# @File    : parallel_sites.py.py
# @annotation    :

from multimer_scan import Mutation
from utils.common import *
from Bio.PDB import PDBParser

def mk_rosetta_resfile(multi_points_mutations: list):
    total_number = len(multi_points_mutations)
    with open('mutfile', 'w+') as mtfile:
        mtfile.write(f'total 1\n{total_number}\n')
        for mutation in multi_points_mutations:
            mtfile.write(" ".join(mutation) + "\n")

THE20 = {'ALA': 0, 'ARG': 1, 'ASN': 2, 'ASP': 3, 'CYS': 4,
         'GLN': 5, 'GLU': 6, 'GLY': 7, 'HIS': 8, 'ILE': 9,
         'LEU': 10, 'LYS': 11, 'MET': 12, 'PHE': 13, 'PRO': 14,
         'SER': 15, 'THR': 16, 'TRP': 17, 'TYR': 18, 'VAL': 19}

pdb = 'test_bak/1NWW.pdb'

parser = PDBParser(QUIET=True)
s = parser.get_structure('s', pdb)
resnum = 1
for i, residue in enumerate(s.get_residues()):
    name_three = residue.get_resname()
    chain = residue.get_parent().id
    resseq = residue.get_id()[1]
    if name_three in long2short:
        resnum += 1
        print(long2short[name_three], chain, resseq, resnum+1)


