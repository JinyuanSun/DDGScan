#!/usr/bin/env python

# https://salilab.org/modeller/manual/node34.html

import os

from modeller import *
import distutils.dir_util


def _3_2_1(x):
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
    try:
        y = d[x]
        return y
    except KeyError:
        return "A"


def getPdbRes(code):
    e = Environ()
    m = Model(e, file=code)
    aln = Alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code + ".seq")
    return code + ".seq"


def generateFillSeq(code, chain, seq=""):

    def outAliFile(seq):
        with open("seqfill.ali", "w+") as seqfill:
            seqfill.write(">P1;seqfill\nsequence:::::::::\n")
            i = 0
            while i + 75 in range(len(seq) + 75):
                if i + 75 <= len(seq):
                    seqfill.write(seq[i: i + 75] + "\n")

                else:
                    seqfill.write(seq[i:] + "*\n")
                i += 75
            seqfill.close()

    def extractSeqFromPDB(code, chain):

        seq = ""
        with open(code + ".pdb", "r") as pdbfile:
            for line in pdbfile:
                if line.startswith("SEQRES"):
                    # print(line)
                    linelist = line.strip().split()
                    if chain == linelist[2]:
                        reslist = linelist[4:]
                        for aa in reslist:
                            seq += _3_2_1(aa)
            pdbfile.close()
        print(f"Extracted Seq: {seq}")
        return seq

    if bool(seq):
        # print(1)
        outAliFile(seq)
        return seq
    else:
        print(
            "[WARNING]: SEQRES information may not exactly match your sequence of you target protein!"
        )
        seq = extractSeqFromPDB(code, chain)
        outAliFile(seq)
        return seq


def align2d(code, chain):
    env = Environ()
    aln = Alignment(env)
    mdl = Model(
        env, file=code, model_segment=("FIRST:%s" % (chain), "LAST:%s" % (chain))
    )
    aln.append_model(
        mdl, align_codes="%s%s" % (code, chain), atom_files="%s.pdb" % (code)
    )
    aln.append(file="seqfill.ali", align_codes="seqfill")
    aln.align2d()
    aln.write(file="pdb-res.ali", alignment_format="PIR")
    aln.write(file="pdb-res.pap", alignment_format="PAP")


def buildModel(code, chain):
    # from modeller import *
    from modeller.automodel import LoopModel, refine  # Load the AutoModel class

    log.none()
    env = Environ()

    # directories for input atom files
    env.io.atom_files_directory = [".", "%s.pdb" % (code)]

    a = LoopModel(
        env, alnfile="pdb-res.ali", knowns="%s%s" % (code, chain), sequence="seqfill"
    )
    a.starting_model = 1
    a.ending_model = 3

    a.loop.starting_model = 1
    a.loop.ending_model = 3
    a.loop.md_level = refine.fast

    a.make()

    ok_models = sorted(
        [x for x in a.loop.outputs if x["failure"] is None],
        key=lambda model: model["molpdf"],
    )
    best = ok_models[0]["name"]

    # print(ok_models)
    # print(best)
    return best


def main(pdb, chain, seq=""):
    # print(bool(seq))
    os.system("cp %s %s.bak" %(pdb, pdb))
    pdb_filename = pdb.split("/")[-1]
    code = pdb_filename.replace(".pdb", "")
    distutils.dir_util.mkpath("modeller")
    os.system("cp %s modeller/" %pdb)
    os.chdir("modeller")
    # code = '4R21'
    # chain = "A"
    seq = generateFillSeq(code, chain, seq)
    align2d(code, chain)
    best = buildModel(code, chain)
    # print(best)
    os.system("cp %s ../%s" % (best, pdb_filename))
    os.chdir("../")
    return seq


if __name__ == "__main__":
    code = "/Users/jsun/RESEARCH/DDGScan/test_bak/1NWW.pdb"
    chain = "A"
    seq = ""
    # if bool(seq):
    #     print("Have sequences")
    # else:
    #     print("None sequence")
    main(code, chain, seq)
