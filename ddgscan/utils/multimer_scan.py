#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/6/17 2:08 AM
# @File    : multimer_scan.py.py
# @annotation    :
from Bio.PDB import PDBParser, NeighborSearch, Atom, Structure, Selection
import pandas as pd
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


class ProteinInterface:
    def __init__(self, pdb: [str, Structure]):
        if type(pdb) == str:
            self.parser = PDBParser(PERMISSIVE=1)
            self.complex = self.parser.get_structure('s', pdb)
        if type(pdb) == Structure:
            self.complex = complex
        self.ca_dict = self.get_ca()
        self.chains = list(self.ca_dict.keys())
        assert len(self.chains) >= 2, f"Only one chain found in {pdb}"
        self.radius = 8
        self.interface_residues = []

    def get_ca(self):
        ca_dict = {}
        for chain in self.complex.get_chains():
            ca_dict[chain.id] = []
            for atom in chain.get_atoms():
                if atom.get_id() == 'CA':
                    ca_dict[chain.id].append(atom)
        return ca_dict

    def get_neighbour_res(self, atom_list: list, center_atom: Atom, radius: float):
        """
        Return residue list of residues' CA within raduis of center atom
        """
        ca_list = []
        for atom in atom_list:
            if atom.get_id() == 'CA':
                ca_list.append(atom)
        nbs = NeighborSearch(ca_list)
        nb_atom_list = nbs.search(center_atom.get_coord(), radius)
        return [atom.get_parent() for atom in nb_atom_list]

    def find_homomultimer_interface(self):
        chain_A_ca = list(self.ca_dict.values())[0]
        other_ca = []
        for ca_atoms in list(self.ca_dict.values())[1:]:
            other_ca += ca_atoms
        #         print(len(other_ca),len(chain_A_ca))
        neighbor_dict = {}

        for atom in chain_A_ca:
            neighbor_dict[atom.get_parent()] = self.get_neighbour_res(other_ca, atom, self.radius)


        for k, v in neighbor_dict.items():
            if len(v) == 0:
                continue
            else:
                self.interface_residues.append(k)
        return neighbor_dict


class Multimerscan:
    # input struc
    # A complementary module for DDGScan, analyze interface positions
    # and make multiple mutations base on beneficial mutations suggested
    def __init__(self, input_pdb: str, seq_file: str, clean: bool = True):
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
        if clean:
            self._remove_nonprotein()

        self.interface = ProteinInterface(self.structure)
        self.interface.find_homomultimer_interface()
        self.mutations = self.read_score_file()

    def _remove_hydrogens(self, structure: Structure):
        # Removes all hydrogens.
        # This code is not suited to work with hydrogens
        for residue in Selection.unfold_entities(structure, 'R'):
            remove = []
            for atom in residue:
                if atom.element == 'H': remove.append(atom.get_id())
                if atom.name == 'OXT': remove.append(atom.get_id())
            for i in remove: residue.detach_child(i)

    def _convert_mse(self, structure: Structure):
        # Changes MSE residues to MET
        for residue in Selection.unfold_entities(structure, 'R'):
            if residue.get_resname() == 'MSE':
                residue.resname = 'MET'
                for atom in residue:
                    if atom.element == 'SE':
                        new_atom = Atom.Atom('SD',
                                             atom.coord,
                                             atom.bfactor,
                                             atom.occupancy,
                                             atom.altloc,
                                             'SD  ',
                                             atom.serial_number,
                                             element='S')
                        residue.add(new_atom)
                        atom_to_remove = atom.get_id()
                        residue.detach_child(atom_to_remove)

    def _remove_water(self, structure: Structure):
        # Removes all water molecules
        residues_to_remove = []
        for residue in Selection.unfold_entities(structure, 'R'):
            if residue.get_resname() == 'HOH':
                residues_to_remove.append(residue)
        for r in residues_to_remove:
            r.get_parent().detach_child(r.get_id())

    def _remove_hetatm(self, structure: Structure):
        # Removes all non-protein molecules
        residues_to_remove = []
        for residue in Selection.unfold_entities(structure, 'R'):
            if residue.id[0] != ' ':
                residues_to_remove.append(residue)
        for r in residues_to_remove:
            r.get_parent().detach_child(r.get_id())

    def _read_pdb(self):
        s = self.parser.get_structure('input', self.input_pdb)
        return s

    def _remove_nonprotein(self):
        # remove water, ligand, solvents, DNA/RNA, and etc.
        self._remove_hydrogens(self.structure)
        self._convert_mse(self.structure)
        self._remove_hetatm(self.structure)

    def read_score_file(self):
        filename = 'Selected_Mutation.csv'
        score_df = pd.read_csv(filename, sep=',', header=0, index_col=0)
        mutations = []
        for mutation in score_df['mutation']:
            mutations.append(Mutation(mutation[0], int(mutation[1:-1]), mutation[-1]))
        return mutations

    def _generate_all_mutations(self, engine: str):
        # make list of mutations

        single_mutations = self.read_score_file()
        interface_residues =  self.interface.interface_residues

        mutations_at_interface = []
        for residue in interface_residues:
            resnum = residue.get_id()[1]
            mutations_at_interface.append(resnum)

        selected_mutations = []
        for mutation in single_mutations:
            if int(mutation.position) in mutations_at_interface:
                selected_mutations.append(mutation)

        if engine.lower() == 'foldx':
            return [mutation.convert2foldx(self.interface.chains[0]) for mutation in selected_mutations]

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
