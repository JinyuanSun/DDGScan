#!/usr/bin/env python
import time


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

    class Model:
        def __init__(self):
            self.chain
            pass



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
        assert AA in d, "%s is not in 20 canonical amino acids!" %(AA)
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

    def read_coords(self):






if __name__ == '__main__':
    pdb_1lve = PDB("1PGA.pdb")
    print(pdb_1lve.read_seqres())
