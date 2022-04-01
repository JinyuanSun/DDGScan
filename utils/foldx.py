#!/usr/bin/env python

# By Jinyuan Sun, Oct, 12, 2021

import distutils.dir_util
import os
import time

import pandas as pd

from utils.common import *


class FoldX:
    def __init__(self, pdbName, path2foldx, numThreads):
        self.pdbname = pdbName
        self.path = path2foldx
        self.threads = numThreads
        self.cutoff: int
        self.result = []

    def repairPDB(self):
        cmd = self.path + "foldx --command=RepairPDB --pdb=" + self.pdbname
        pdbname = self.pdbname
        # print(cmd)
        FoldX_out = os.popen(cmd).read()
        with open(".foldx_repair.log", "w+") as outfile:
            outfile.write(FoldX_out)
            outfile.close()
        return pdbname.replace(".pdb", "_Repair.pdb")

    def calScore(self, wild, resNum, mutation, pdbfile, jobID):
        fxout_name = jobID + "/Dif_" + pdbfile.replace(".pdb", ".fxout")
        # print(fxout_name)
        df = pd.read_table(fxout_name, sep="\t", skiprows=8)
        score = round(df["total energy"].mean(), 4)
        sd = round(df["total energy"].std(), 4)
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

        with open("individual_list.txt", "w+") as indFile:
            indFile.write(wild + chain + str(resNum) + mutation + ";")
            indFile.close()
        cmd1 = "cp ../../" + pdbfile + " ./"
        os.popen(cmd1)
        cmd2 = (
                self.path
                + "foldx --command=BuildModel --numberOfRuns="
                + numOfRuns
                + " --mutant-file=individual_list.txt --pdb="
                + pdbfile
                + " 1>/dev/null"
        )

        starttime = time.time()
        os.system(cmd2)
        finishtime = time.time()
        print(
            "[DEBUG]: FoldX mutation %s_%s_%s took %f seconds."
            % (wild, resNum, mutation, finishtime - starttime)
        )

        os.chdir("../../")
        # return self.calScore(self, pdbfile)  # need self?


class foldx_binder:
    def __init__(self):
        pass

    @staticmethod
    def repair_pdb(pdb_file):
        cmd = "foldx --command=RepairPDB --pdb=" + pdb_file
        repair_out = os.popen(cmd).read()
        with open(".foldx_repair.log", "w+") as outfile:
            outfile.write(repair_out)
            outfile.close()
        return pdb_file.replace(".pdb", "_Repair.pdb")

    @staticmethod
    def cal_score(wild, resNum, mutation, pdbfile):
        fxout_name = "Dif_" + pdbfile.replace(".pdb", ".fxout")
        # print(fxout_name)
        df = pd.read_table(fxout_name, sep="\t", skiprows=8)
        score = round(df["total energy"].mean(), 4)
        min_score = round(df["total energy"].min(), 4)
        sd = round(df["total energy"].std(), 4)
        # result.append(["_".join([wild, str(resNum), mutation]), score, sd])
        return ["_".join([wild, str(resNum), mutation]), str(score), str(min_score), str(sd)]

    @staticmethod
    def run_one_job(varlist: list):
        pdb_file, wild, chain, mutation, position, job_id, numOfRuns = varlist
        # mutation_name = "_".join([wild, str(resNum), mutation])
        path_job_id = FOLDX_JOBS_DIR + job_id
        distutils.dir_util.mkpath(path_job_id)
        os.chdir(path_job_id)
        with open("individual_list.txt", "w+") as indFile:
            indFile.write(wild + chain + str(position) + mutation + ";")
            indFile.close()
        cmd1 = "cp ../../" + pdb_file + " ./"
        os.popen(cmd1)
        cmd2 = (
                "foldx --command=BuildModel --numberOfRuns="
                + str(numOfRuns)
                + " --mutant-file=individual_list.txt --pdb="
                + pdb_file
                + " 1>/dev/null"
        )
        starttime = time.time()
        os.system(cmd2)
        results = foldx_binder.cal_score(wild, position, mutation, pdb_file)
        cp_files(job_id, pdb_file, numOfRuns)
        finishtime = time.time()
        print(
            "FoldX mutation %s_%s_%s took %f seconds."
            % (wild, position, mutation, finishtime - starttime)
        )

        os.chdir("../../")
        return results


def cp_files(job_id, pdb, num_of_runs):
    pdb_id = pdb.replace(".pdb", "")
    distutils.dir_util.mkpath('../inspection')
    names = [pdb]
    os.system(f"cp {pdb} ../inspection")
    for i in range(int(num_of_runs)):
        fxout_pdb_name = pdb_id + "_1_" + str(i) + ".pdb"
        os.system(f"cp {fxout_pdb_name} ../inspection/{pdb_id}_{job_id}_{i}.pdb")
        names.append(f'{pdb_id}_{job_id}_{i}.pdb')
    with open(f'../inspection/{job_id}.pml', 'w+') as pml:
        for name in names:
            pml.write(f'load {name}\n')
        pml.write(f'select mutation, resi {job_id.split("_")[1]}\n')
        pml.write("select br. all within 6 of mutation\n")
        pml.write("show sticks, sele\n")
        pml.write('color indium, mutation AND elem C AND sc.')
        pml.close()
