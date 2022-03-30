#!/usr/bin/env python
from .io import Protein


def judge(userSeq, seq, resNumList):
    listLen = int(resNumList[-1]) - int(resNumList[0]) + 1
    resLen = len(seq)
    # print(listLen, resLen)
    if listLen != resLen:
        return userSeq
    else:
        return 0  # no ncAA and gap found


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
