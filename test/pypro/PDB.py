#!/usr/bin/env python
import time
import structure


class PDB:
    def __init__(self, filename=''):
        import time
        # self.Annotation = PDB.Annotation()
        self.Model = PDB.Model
        self.pdbfilelist = []
        self.pdbdict = {}
        if filename:
            self.pdbfilelist = self.read_in_file(filename)
            # self.SEQRES = self.read_seqres()


    def timing(self):
        pass

    def three2oneAA(self, AA):
        d = {
            "CYS": "C",
            "ASP": "D",
            "SER": "S",
            "GLN": "Q",
            "LYS": "K",
            "ILE": "I",
            "PRO": "P",
            "THR": "T",
            "PHE": "F",
            "ASN": "N",
            "GLY": "G",
            "HIS": "H",
            "LEU": "L",
            "ARG": "R",
            "TRP": "W",
            "ALA": "A",
            "VAL": "V",
            "GLU": "E",
            "TYR": "Y",
            "MET": "M",
        }
        assert AA in d, "%s is not in 20 canonical amino acids!" % (AA)
        return d[AA]

    def read_in_file(self, filename):
        with open(filename, 'r') as pdbfile:
            for line in pdbfile:
                self.pdbfilelist.append(line.replace("\n", ""))
            pdbfile.close()
        return self.pdbfilelist

    def read_seqres(self):
        SEQRES_dict = {}
        for line in self.pdbfilelist:
            if line[:6] == "SEQRES":
                chain = line[11]
                reslist = line[19:].strip().split()
                subseq = ""
                for AA in reslist:
                    subseq += self.three2oneAA(AA)
                if chain in SEQRES_dict.keys():
                    SEQRES_dict[chain] += subseq
                else:
                    SEQRES_dict[chain] = subseq
        return SEQRES_dict

    def read_atoms(self):
        # ATOM      1  N   MET A   1      26.778  34.213  35.880  1.00 14.61           N
        # ATOM      2  CA  MET A   1      26.659  32.769  36.242  1.00 16.66           C
        # ATOM      3  C   MET A   1      27.468  31.927  35.268  1.00 16.16           C
        # ATOM      4  O   MET A   1      27.699  32.342  34.110  1.00 15.79           O
        for line in self.pdbfilelist:
            resnum_list = []

            res_dict = {}
            if line[:4] == "ATOM":
                atom = aminoacid.Atom(line)
                if atom.resSeq in resnum_list:
                    res_dict[atom.resSeq][atom.name] = atom
                else:
                    res_dict[atom.resSeq] = {atom.name: atom}
            return







if __name__ == '__main__':
    pdb_1lve = PDB("1PGA.pdb")
    print(pdb_1lve.read_seqres())
