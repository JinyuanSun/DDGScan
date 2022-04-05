#!/usr/bin/env python
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s   %(levelname)s   %(message)s')


class Protein:
    def __init__(self, pdbname, chain):
        self.pdbname = pdbname
        self.chain = chain
        self.seq = ''
        self.resNumList = []

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
        y = d.get(x)
        assert y, f"{x} dose not belong to 20 canonical animo acids!"
        return y

    def pdb2seq(self):
        with open(self.pdbname) as pdbfile:
            for line in pdbfile:
                if "ATOM" == line[0:6].replace(" ", ""):
                    if self.chain == line[21].replace(" ", ""):
                        if line[12:16].replace(" ", "") == "CA":
                            if line[16] == "B":
                                continue
                            else:
                                self.seq += self._3_2_1(line[17:20].replace(" ", ""))
                                self.resNumList.append(
                                    int(line[22:26].replace(" ", ""))
                                )
        pdbfile.close()
        return self.seq, self.resNumList


def judge(userSeq, seq, resNumList):
    # print(resNumList, resNumList)
    listLen = int(resNumList[-1]) - int(resNumList[0]) + 1
    resLen = len(seq)
    # print(listLen, resLen)
    if userSeq:
        if listLen != resLen:
            return userSeq
        else:
            if listLen != len(userSeq):
                print("Missing at N- or C-terminal is detected, however grape will not build it!")
            return 0
    else:  # no sequence provided
        if listLen != resLen:
            return "chain break found"
        else:
            return 0  # no ncAA and gap found


def main(pdb, chain, userSeq):
    seq, resNumList = Protein(pdb, chain).pdb2seq()
    judge_result = judge(userSeq, seq, resNumList)
    if judge_result != 0:
        logging.error("Chain break detected! Transfer to AlphaFold.")
        exit()


if __name__ == "__main__":
    import sys

    pdb = sys.argv[1]
    chain = sys.argv[2]
    userSeq = sys.argv[3]
    main(pdb, chain, userSeq)
