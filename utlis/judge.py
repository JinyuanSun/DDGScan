#!/usr/bin/env python


class Protein:
    def __init__(self, pdbname, chain):
        self.pdbname = pdbname
        self.chain = chain
        self.seq, self.resNumList = self.pdb2seq()

    def _3_2_1(self, x):
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
        y = d[x]
        return y

    def pdb2seq(self):
        seq = ""
        resNumList = []
        with open(self.pdbname) as pdbfile:
            for line in pdbfile:
                if "ATOM" == line[0:6].replace(" ", ""):
                    if self.chain == line[21].replace(" ", ""):
                        if line[12:16].replace(" ", "") == "CA":
                            if line[16] == "B":
                                # print(line)
                                continue
                            else:
                                seq += self._3_2_1(line[17:20].replace(" ", ""))
                                resNumList.append(int(line[22:26].replace(" ", "")))
        pdbfile.close()
        return seq, resNumList


def judge(userSeq, seq, resNumList):
    listLen = int(resNumList[-1]) - int(resNumList[0]) + 1
    resLen = len(seq)
    # print(listLen, resLen)
    if listLen != resLen:
        return userSeq
    else:
        return 0

def main(pdb, chain, userSeq):
    userPdb = Protein(pdb, chain)
    return judge(userSeq, userPdb.seq, userPdb.resNumList)



if __name__ == "__main__":
    import sys

    pdb = sys.argv[1]
    chain = sys.argv[2]
    userSeq = sys.argv[3]
    userPdb = Protein(pdb, chain)
    judge_result = judge(userSeq, userPdb.seq, userPdb.resNumList)
    # print(userSeq, userPdb.seq, userPdb.resNumList)
    if judge_result != 0:
        print("Chain break detected! Transfer to AlphaFold.")
