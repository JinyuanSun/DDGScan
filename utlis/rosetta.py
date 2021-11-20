#!/usr/bin/env python

# By Jinyuan Sun, Oct, 13, 2021

import os
import pandas as pd
import numpy as np
import subprocess
import time

class Rosetta:
    def __init__(self, pdbName, relax_num, numThreads, exe='cartesian_ddg.static.linuxgccrelease', rosettadb= "/opt/rosetta_bin_linux_2021.16.61629_bundle/main/database"):
        self.exe = exe
        self.pdbname = pdbName
        self.relax_num = relax_num
        self.threads = numThreads
        # self.relaxedpdb = pdbName #for test
        self.rosettadb = rosettadb
        self.relaxedpdb: str  # for test

        #self.cutoff = cutOff
        self.result = []

    def relax(self):
        try:
            os.mkdir('rosetta_relax')
            os.system("cp  " + self.pdbname + " rosetta_relax/")
            os.chdir('rosetta_relax')
        except FileExistsError:
            os.system("cp  " + self.pdbname + " rosetta_relax/")
            os.chdir('rosetta_relax')
            pass

        with open("cart2.script","w+") as cart2:
            cart2.write("switch:cartesian\n")
            cart2.write("repeat 2\n")
            cart2.write("ramp_repack_min 0.02  0.01     1.0  50\n")
            cart2.write("ramp_repack_min 0.250 0.01     0.5  50\n")
            cart2.write("ramp_repack_min 0.550 0.01     0.0 100\n")
            cart2.write("ramp_repack_min 1     0.00001  0.0 200\n")
            cart2.write("accept_to_best\n")
            cart2.write("endrepeat")
            cart2.close()
        relax_cmd = "".join(
            ["mpirun -n " + str(self.threads) + " relax.mpi.linuxgccrelease -s " + self.pdbname + " -use_input_sc",
            " -constrain_relax_to_start_coords -ignore_unrecognized_res",
            " -nstruct " + str(self.relax_num),
            " -relax:coord_constrain_sidechains",
            " -relax:cartesian -score:weights ref2015_cart ",
            " -relax:min_type lbfgs_armijo_nonmonotone",
            " -relax:script cart2.script 1>/dev/null && sort -nk2 score.sc |head -n 1|awk '{print$22}'"]
        )
        print("==" * 20)
        print(" Relaxing your Protein: ")
        # os.system(relax_cmd)
        relaxed_pdb_name = os.popen(relax_cmd).read()
        print(" Finished relax! ")
        print("==" * 20)
        relaxed_pdb_name = os.popen("sort -nk2 score.sc |head -n 1|awk '{print$22}'").read()
        self.relaxedpdb = relaxed_pdb_name.replace("\n", "") + ".pdb"
        os.chdir("../")
        return relaxed_pdb_name.replace("\n", "") + ".pdb"

    def read_rosetta_ddgout(self, rosettaddgfilename):
        ddg_dict = {}

        ddg_array = []
        with open(rosettaddgfilename) as rosettaddg:
            for line in rosettaddg:
                #print(line.split(":")[2])
                if line.split(":")[2].strip() == "WT":
                    dg_ref = float(line.split(":")[3][1:10])
                else:
                    ddg = float(line.split(":")[3][1:10]) - dg_ref
                    ddg_array.append(ddg)
            rosettaddg.close()

        return round(np.array(ddg_array).mean(),4), round(np.array(ddg_array).std(),4)

    def runOneJob(self, varlist: list):
        wild, mutation, resNum, jobID = varlist
        try:
            os.mkdir(jobID)
            os.chdir(jobID)
        except FileExistsError:
            os.chdir(jobID)

        # os.popen('cp ../../rosetta_relax/' + self.relaxedpdb + ' ./')
        os.system('cp ../../rosetta_relax/' + self.relaxedpdb + ' ./')
        with open('mtfile', 'w+') as mtfile:
            mtfile.write("total 1\n")
            mtfile.write("1\n")
            mtfile.write(wild + " " + str(resNum) + " " + mutation)
            mtfile.close()

        argument_list = [self.exe, "-database", self.rosettadb,
                         "-use_input_sc",
                         "-s", self.relaxedpdb,
                         "-ddg:mut_file", "mtfile",
                         "-ddg:iterations", "3",
                         "-ddg::cartesian",
                         "-ddg::dump_pdbs", "true",
                         "-ddg:bbnbrs", "1",
                         "-score:weights", "ref2015_cart",
                         "-relax:cartesian",
                         "-relax:min_type", "lbfgs_armijo_nonmonotone",
                         "-flip_HNQ",
                         "-crystal_refine",
                         "-fa_max_dis",
                         "9.0",
                         "1>/dev/null"]
        cartddg_cmd = " ".join(argument_list)
        os.system(cartddg_cmd)
        print(cartddg_cmd)
        os.chdir("../../")
        # return pid, '_'.join([wild, str(trueResNum), mutation])


if __name__ == '__main__':
    pdbname = '1PGA.pdb'
    chain = 'A'
    threads = 24
    relax_num = 200
    prot = Rosetta(pdbname, relax_num, 'numThreads')
    score, std = prot.read_rosetta_ddgout('rosetta_jobs/0_1_/mtfile.ddg')
    print(score, std)
    #relax_pdb_name = prot.relax(pdbname, threads)
    #print("Using: " + relax_pdb_name)
    #os.mkdir('rosetta_jobs')
    #prot.runOneJob(relax_pdb_name, "M", chain, "Q", 1, 'rosetta_jobs/0_1/')
