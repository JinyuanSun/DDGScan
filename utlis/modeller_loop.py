#!/usr/bin/env python

# https://salilab.org/modeller/manual/node34.html

import os

from modeller import *


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
                    linelist = line.strip().split()
                    if chain == linelist[2]:
                        reslist = linelist[4:]
                        for aa in reslist:
                            seq += _3_2_1(aa)
            pdbfile.close()
        return seq

    if bool(seq):
        outAliFile(seq)
    else:
        print(
            "[WARNING]: SEQRES information may not exactly match your sequence of you target protein!"
        )
        seq = extractSeqFromPDB(code, chain)
        outAliFile(seq)


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
    code = pdb.replace(".pdb", "")
    os.mkdir("modeller")
    os.chdir("modeller")
    os.system("cp ../%s ./" % (pdb))
    # code = '4R21'
    # chain = "A"
    generateFillSeq(code, chain, seq)
    align2d(code, chain)
    best = buildModel(code, chain)
    os.system("cp %s ../%s_modfixed.pdb" % (best, code))
    os.chdir("../")
    return "%s_modfixed.pdb" % (code)


if __name__ == "__main__":
    code = "4R21"
    chain = "A"
    seq = ""
    main(code, chain, seq)
