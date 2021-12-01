#!/usr/bin/env python

import pandas as pd
import os
import utlis.foldx as foldx
import utlis.rosetta as rosetta

import utlis.io as io
from joblib import Parallel, delayed

pdb = '6LR7.pdb'

prot_rosetta = rosetta.Rosetta(pdb, 50, 4)
# relaxed_prot = prot_rosetta.relax()

prot = io.Protein(pdb, 'A')
seq, resNumList = io.Protein.pdb2seq(prot)

try:
    os.mkdir('rosetta_jobs')
except FileExistsError:
    pass
# all_results = []
job_list = []
for i, res in enumerate(seq):
    # i = int(i)
    resNum = resNumList[i]
    wild = res
    for j, aa in enumerate('QWERTYIPASDFGHKLCVNM'):
        if aa != wild:
            jobID = "rosetta_jobs/" + "_".join([wild, str(resNum), aa])
            job_list.append([wild, aa, str(i + 1), jobID])


Parallel(n_jobs=4)(delayed(print)(var) for var in job_list)

