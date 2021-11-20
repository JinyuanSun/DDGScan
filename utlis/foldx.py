#!/usr/bin/env python

# By Jinyuan Sun, Oct, 12, 2021

import os
import pandas as pd


class FoldX:
    def __init__(self, pdbName, path2foldx, numThreads, cutOff):
        self.pdbname = pdbName
        self.path = path2foldx
        self.threads = numThreads
        self.cutoff = cutOff
        self.result = []

    def repairPDB(self):
        cmd = self.path + 'foldx --command=RepairPDB --pdb=' + self.pdbname
        pdbname = self.pdbname
        # print(cmd)
        print('==' * 20)
        print('Stage 1: Repair PDB file using FoldX')
        print('==' * 20)
        print('[INFO]: Running FoldX RepairPDB.')

        FoldX_out = os.popen(cmd).read()
        with open('.foldx_repair.log', 'w+') as outfile:
            outfile.write(FoldX_out)
            outfile.close()

        print('[INFO]: FoldX RepairPDB Done!')
        return pdbname.replace(".pdb", '_Repair.pdb')

    def calScore(self, wild, resNum, mutation, pdbfile, jobID):
        fxout_name = jobID + "/Dif_" + pdbfile.replace(".pdb", '.fxout')
        # print(fxout_name)
        df = pd.read_table(fxout_name, sep='\t', skiprows=8)
        score = round(df['total energy'].mean(), 4)
        sd = round(df['total energy'].std(), 4)
        self.result.append(["_".join([wild, str(resNum), mutation]), score, sd])
        return ["_".join([wild, str(resNum), mutation]), score, sd]

    def runOneJob(self, varlist: list):
        pdbfile, wild, chain, mutation, resNum, jobID, numOfRuns = varlist
        # mutation_name = "_".join([wild, str(resNum), mutation])
        try:
            os.mkdir(jobID)
            os.chdir(jobID)
        except FileExistsError:
            os.chdir(jobID)

        with open('individual_list.txt', 'w+') as indFile:
            indFile.write(wild + chain + str(resNum) + mutation + ";")
            indFile.close()
        cmd1 = "cp ../../" + pdbfile + " ./"
        os.popen(cmd1)
        cmd2 = self.path + "foldx --command=BuildModel --numberOfRuns=" + numOfRuns + \
               " --mutant-file=individual_list.txt --pdb=" + pdbfile + " 1>/dev/null"
        os.system(cmd2)
        os.chdir("../../")
        # return self.calScore(self, pdbfile)  # need self?
