#!/usr/bin/env python

# By Jinyuan Sun, Oct, 12, 2021
import argparse
class Protein:
    def __init__(self, pdbname, chain):
        self.pdbname = pdbname
        self.chain = chain
        self.seq = ''
        self.resNumList = []

    def _3_2_1(self, x):
        d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
        y = d[x]
        return y

    def pdb2seq(self):
        with open(self.pdbname) as pdbfile:
            for line in pdbfile:
                if 'ATOM' == line[0:6].replace(" ",""):
                    if self.chain == line[21].replace(" ",""):
                        if line[12:16].replace(" ","") == "CA":
                            self.seq += self._3_2_1(line[17:20].replace(" ",""))
                            self.resNumList.append(int(line[22:26].replace(" ","")))
        pdbfile.close()
        return self.seq, self.resNumList

class Parser:

    def __init__(self):
        pass
    def get_args():
        parser = argparse.ArgumentParser(description='Run GRAPE')
        parser.add_argument("-pdb", '--pdb', help="Input PDB")
        parser.add_argument("-chain", '--chain', help="Input PDB Chain to do in silico DMS", default="A")
        parser.add_argument("-cpu", '--threads', help="Number of threads to run FoldX, Rosetta and HHblits", default=16)

        # parser.add_argument("-r", '--ratio', help="Select by ratio mode",default=False)

        parser.add_argument("-fc", '--foldx_cutoff',help="Cutoff of foldx ddg", default=1.5)

        parser.add_argument("-mode", '--mode',help="Run, Rerun or analysis")


        args = parser.parse_args()

        return args


if __name__ == '__main__':
    pdbname = '6JTT.pdb'
    chain = 'A'
    prot = Protein(pdbname, chain)
    seq, resNumList = Protein.pdb2seq(prot)
    print(seq, resNumList)


