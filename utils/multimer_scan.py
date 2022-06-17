#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/6/17 2:08 AM
# @File    : multimer_scan.py.py
# @annotation    :
from Bio.PDB import PDBParser
from common import *

class Mutation:
    def __init__(self, w: str, p: int, m: str):
        self.wildtype = w
        self.position = p
        self.mutation = m

    def convert2rosetta(self, abs_number: int):
        return [self.wildtype, str(abs_number), self.mutation]

    def convert2foldx(self, chain: str):
        return [self.wildtype, chain, str(self.position), self.mutation]


class Multimerscan:
    def __init__(self, input_pdb, seq_file):
        self.input_pdb = input_pdb
        self.input_seq = seq_file
        self.parser = PDBParser(PERMISSIVE=True)
        self.structure = self._read_pdb()
        self.structure2chain_dict()
        if self.input_seq:
            self._read_fasta(self.input_seq)
            self.homo_chains = self.detect_homo(self.fasta_chain_dict)
        else:
            self.homo_chains = self.detect_homo(self.pdb_chain_dict)

    def _read_pdb(self):
        s = self.parser.get_structure('input', self.input_pdb)
        return s

    def _remove_nonprotein(self):
        # remove water, ligand, solvents, DNA/RNA, and etc.
        pass

    def _position_mapping(self):
        # map pdb numbering to reference sequence
        pass

    def _generate_all_mutations(self):
        # make list of mutations
        pass

    def _read_fasta(self, seq_file):
        self.fasta_chain_dict = {}
        with open(seq_file, 'r') as fasta:
            for line in fasta:
                if line.startswith(">"):
                    chain_name = line.replace(">", "").replace("\n", "")
                    self.fasta_chain_dict[chain_name] = ''
                else:
                    self.fasta_chain_dict[chain_name] += line.replace("\n", "")

    def structure2chain_dict(self):
        self.pdb_chain_dict = {}
        for i, residue in enumerate(self.structure.get_residues()):
            resname = residue.get_resname()
            if resname in long2short:
                chain = residue.get_parent().id
                if chain in self.pdb_chain_dict:
                    self.pdb_chain_dict[chain] += long2short(resname)
                else:
                    self.pdb_chain_dict[chain] = long2short(resname)

    def detect_homo(self, chain_dict):
        homo_chains = {}
        for i, (chain_1, seq_1) in enumerate(chain_dict.items()):
            picked = "".join(["".join(x) for x in homo_chains.values()])
            if chain_1 in picked:
                pass
            else:
                homo_chains[chain_1] = []
                for j, (chain_2, seq_2) in enumerate(chain_dict.items()):
                    if j > i and seq_1 == seq_2:
                        homo_chains[chain_1].append(chain_2)
        return homo_chains

    def generate_mutations(self, chain_dict: dict, homo_chains: dict):
        pass




