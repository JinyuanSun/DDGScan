#!/usr/bin/env python

# By Jinyuan Sun, Oct, 13, 2021

import distutils.dir_util
import os
import time
from shutil import which
import numpy as np
import utils.common as common
# from .common import *


class Rosetta:
    def __init__(
            self,
            pdbName,
            relax_num,
            numThreads,
            exe,
            rosettadb,
    ):
        self.exe = exe
        self.pdbname = pdbName
        self.relax_num = relax_num
        self.threads = numThreads
        # self.relaxedpdb = pdbName #for test
        self.rosettadb = rosettadb
        self.relaxedpdb: str  # for test
        self.result = []

    def relax(self):
        distutils.dir_util.mkpath(common.ROSETTA_RELAX_DIR)
        os.system("cp  " + self.pdbname + " " + common.ROSETTA_RELAX_DIR)
        os.chdir(common.ROSETTA_RELAX_DIR)

        with open("cart2.script", "w+") as cart2:
            cart2.write("switch:cartesian\n")
            cart2.write("repeat 2\n")
            cart2.write("ramp_repack_min 0.02  0.01     1.0  50\n")
            cart2.write("ramp_repack_min 0.250 0.01     0.5  50\n")
            cart2.write("ramp_repack_min 0.550 0.01     0.0 100\n")
            cart2.write("ramp_repack_min 1     0.00001  0.0 200\n")
            cart2.write("accept_to_best\n")
            cart2.write("endrepeat")
            cart2.close()

        relax_threads = min([int(self.threads), int(self.relax_num)])
        relax_cmd = "".join(
                [
                    "mpirun --allow-run-as-root -n "
                    + str(relax_threads)
                    + " relax.mpi.linuxgccrelease -s "
                    + self.pdbname
                    + " -use_input_sc",
                    " -constrain_relax_to_start_coords -ignore_unrecognized_res",
                    " -nstruct " + str(self.relax_num),
                    " -relax:coord_constrain_sidechains",
                    " -relax:cartesian -score:weights ref2015_cart ",
                    " -relax:min_type lbfgs_armijo_nonmonotone",
                    " -relax:script cart2.script 1>/dev/null && sort -nk2 score.sc |head -n 1|awk '{print$22}'",
                ]
            )
        if which('relax.mpi.linuxgccrelease'):
            print("==" * 20)
            print(" Relaxing your Protein: ")
            # os.system(relax_cmd)
            relaxed_pdb_name = os.popen(relax_cmd).read()
            print(" Finished relax! ")
            print("==" * 20)
            relaxed_pdb_name = os.popen(
                "sort -nk2 score.sc |head -n 1|awk '{print$22}'"
            ).read()
            self.relaxedpdb = relaxed_pdb_name.replace("\n", "") + ".pdb"
        else:
            print("==" * 20)
            print("Relax skipped!")
            print("==" * 20)
            relaxed_pdb_name = self.pdbname.replace('.pdb', "")
            self.relaxedpdb = relaxed_pdb_name.replace("\n", "") + ".pdb"
        os.chdir("../")
        return relaxed_pdb_name.replace("\n", "") + ".pdb"

    def fast_relax(self):
        distutils.dir_util.mkpath(common.ROSETTA_RELAX_DIR)
        os.system("cp  " + self.pdbname + " " + common.ROSETTA_RELAX_DIR)
        os.chdir(common.ROSETTA_RELAX_DIR)
        relax_cmd = f"relax.mpi.linuxgccrelease -s {self.pdbname} 1>/dev/null"
        x = os.popen(relax_cmd).read()
        self.relaxedpdb = self.pdbname.replace(".pdb", "_0001.pdb")
        os.chdir("../")
        return self.pdbname.replace(".pdb", "_0001.pdb")

    def read_rosetta_ddgout(self, rosettaddgfilename):

        ddg_array = []
        with open(rosettaddgfilename) as rosettaddg:
            for line in rosettaddg:
                if line.split(":")[2].strip() == "WT":
                    dg_ref = float(line.split(":")[3][1:10])
                else:
                    ddg = float(line.split(":")[3][1:10]) - dg_ref
                    ddg_array.append(ddg)
            rosettaddg.close()

        return [
            "{:.4f}".format(np.array(ddg_array).mean()),
            "{:.4f}".format(np.array(ddg_array).std())
        ]

    def read_ddg_monomer_out(self, ddg_monomer_file):
        ddg = float(open(ddg_monomer_file, 'r').readlines()[1].split()[2])
        return ["{:.4f}".format(ddg), "{:.4f}".format(ddg)]

    def runOneJob(self, varlist: list):
        wild, mutation, resNum, jobID = varlist
        distutils.dir_util.mkpath(jobID)
        os.chdir(jobID)

        os.system("cp ../../" + common.ROSETTA_RELAX_DIR +
                  self.relaxedpdb + " ./")
        with open("mtfile", "w+") as mtfile:
            mtfile.write("total 1\n")
            mtfile.write("1\n")
            mtfile.write(wild + " " + str(resNum) + " " + mutation + "\n")
            mtfile.close()

        argument_list = [
            self.exe,
            "-database",
            self.rosettadb,
            "-use_input_sc",
            "-s",
            self.relaxedpdb,
            "-ddg:mut_file",
            "mtfile",
            "-ddg:iterations",
            "3",
            "-ddg::cartesian",
            "-ddg::dump_pdbs",
            "true",
            "-ddg:bbnbrs",
            "1",
            "-score:weights",
            "ref2015_cart",
            "-relax:cartesian",
            "-relax:min_type",
            "lbfgs_armijo_nonmonotone",
            "-flip_HNQ",
            "-crystal_refine",
            "-fa_max_dis",
            "9.0",
            "1>/dev/null",
        ]
        cartddg_cmd = " ".join(argument_list)

        starttime = time.time()
        os.system(cartddg_cmd)
        finishtime = time.time()
        print(
            "[DEBUG]: Rosetta mutation %s_%s_%s took %f seconds."
            % (wild, resNum, mutation, finishtime - starttime)
        )
        # print(cartddg_cmd)
        os.chdir("../../")

    def pmut_scan(self, relaxed_pdb):
        if os.path.isfile("pmut.out"):
            pass
        else:
            pmut_scan_exe = (
                os.popen("which pmut_scan_parallel.mpi.linuxgccrelease")
                .read()
                .replace("\n", "")
            )
            rosettadb = "/".join(pmut_scan_exe.split("/")[:-3]) + "/database/"
            arg_list = [
                "mpirun",
                "-np",
                str(self.threads),
                pmut_scan_exe,
                "-database",
                rosettadb,
                "-s",
                relaxed_pdb,
                "-ex1",
                "-ex2",
                "-extrachi_cutoff 1",
                "-use_input_sc",
                "-ignore_unrecognized_res",
                "-no_his_his_pairE",
                "-multi_cool_annealer",
                "10",
                "-mute",
                "basic",
                "core" ">",
                "pmut.out && ls pmut.out",
            ]
            print("[INFO]: Running pmut_scan_parallel")
            pmut_start = time.time()
            print(" ")
            os.system(" ".join(arg_list))
            pmut_end = time.time()
            pmut_time = pmut_end - pmut_start
            return pmut_time

    def pmut_scan_analysis(self, pmutoutfile):
        with open(pmutoutfile) as pmut_out:

            i = -1
            start_line = 0
            for line in pmut_out:
                i += 1
                line = line.replace("\x1b[0m", "")
                if line.endswith(
                        "mutation   mutation_PDB_numbering   average_ddG   average_total_energy\n"
                ):
                    start_line = 1
                if "protocol took" in line:
                    start_line = 0
                if start_line == 1:
                    mutinfo = line.replace("\n", "").split(")")[1].split()
                    if mutinfo[2] == "average_ddG":
                        with open(common.ROSETTA_SCORE_FILE, "w") as scorefile:
                            scorefile.write(
                                "#Score file formatted by GRAPE from Rosetta.\n#mutation\tscore\tstd\n"
                            )
                            scorefile.close()
                    else:
                        with open(common.ROSETTA_SCORE_FILE, "a+") as scorefile:
                            mut = mutinfo[0].split("-")[1]
                            scorefile.write(
                                "_".join([mut[0], mut[1:-1], mut[-1]])
                                + "\t"
                                + mutinfo[2]
                                + "\t"
                                + "0\n"
                            )
                            scorefile.close()
        pmut_out.close()


class rosetta_binder:
    def __init__(self):
        pass

    @staticmethod
    def relax(pdbname, threads, relax_num):
        distutils.dir_util.mkpath(common.ROSETTA_RELAX_DIR)
        os.system("cp  " + pdbname + " " + common.ROSETTA_RELAX_DIR)
        os.chdir(common.ROSETTA_RELAX_DIR)

        with open("cart2.script", "w+") as cart2:
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
            [
                "mpirun -n "
                + str(threads)
                + " relax.mpi.linuxgccrelease -s "
                + pdbname
                + " -use_input_sc",
                " -constrain_relax_to_start_coords -ignore_unrecognized_res",
                " -nstruct " + str(relax_num),
                " -relax:coord_constrain_sidechains",
                " -relax:cartesian -score:weights ref2015_cart ",
                " -relax:min_type lbfgs_armijo_nonmonotone",
                " -relax:script cart2.script 1>/dev/null && sort -nk2 score.sc |head -n 1|awk '{print$22}'",
            ]
        )
        print("==" * 20)
        print(" Relaxing your Protein: ")
        # os.system(relax_cmd)
        relaxed_pdb_name = os.popen(relax_cmd).read()
        print(" Finished relax! ")
        print("==" * 20)
        relaxed_pdb_name = os.popen(
            "sort -nk2 score.sc |head -n 1|awk '{print$22}'"
        ).read()
        relaxedpdb = relaxed_pdb_name.replace("\n", "") + ".pdb"
        os.chdir("../")
        return relaxedpdb

    @staticmethod
    def read_rosetta_ddgout(rosettaddgfilename, wild, mutation, resNum):
        ddg_array = []
        with open(rosettaddgfilename, 'r') as rosettaddg:
            for line in rosettaddg:
                # print(line.split(":")[2])
                if line.split(":")[2].strip() == "WT":
                    dg_ref = float(line.split(":")[3][1:10])
                else:
                    ddg = float(line.split(":")[3][1:10]) - dg_ref
                    ddg_array.append(ddg)
            rosettaddg.close()

        return [
            "_".join([wild, str(resNum), mutation]),
            str(round(np.array(ddg_array).mean(), 4)),
            str(round(min(np.array(ddg_array)), 4)),
            str(round(np.array(ddg_array).std(), 4))
        ]

    @staticmethod
    def read_ddg_monomer_out(ddg_monomer_file, wild, mutation, resNum):
        ddg = float(open(ddg_monomer_file, 'r').readlines()[1].split()[2])
        return ["_".join([wild, str(resNum), mutation]),
                "{:.4f}".format(ddg),
                "{:.4f}".format(ddg),
                "0.000"]

    @staticmethod
    def run_one_job(varlist: list):
        wild, mutation, resNum, jobID, relaxedpdb, exe, rosettadb = varlist
        path_job_id = common.ROSETTA_JOBS_DIR + jobID
        distutils.dir_util.mkpath(path_job_id)
        os.chdir(path_job_id)

        os.system("cp ../../" + common.ROSETTA_RELAX_DIR + relaxedpdb + " ./")
        with open("mtfile", "w+") as mtfile:
            mtfile.write("total 1\n")
            mtfile.write("1\n")
            mtfile.write(wild + " " + str(resNum) + " " + mutation + "\n")
            mtfile.close()

        argument_list = [
            exe,
            "-database",
            rosettadb,
            "-use_input_sc",
            "-s",
            relaxedpdb,
            "-ddg:mut_file",
            "mtfile",
            "-ddg:iterations",
            "3",
            "-ddg::cartesian",
            "-ddg::dump_pdbs",
            "true",
            "-ddg:bbnbrs",
            "1",
            "-score:weights",
            "ref2015_cart",
            "-relax:cartesian",
            "-relax:min_type",
            "lbfgs_armijo_nonmonotone",
            "-flip_HNQ",
            "-crystal_refine",
            "-fa_max_dis",
            "9.0",
            "1>/dev/null",
        ]
        cartddg_cmd = " ".join(argument_list)

        starttime = time.time()
        os.system(cartddg_cmd)
        finishtime = time.time()
        print(
            "[DEBUG]: Rosetta mutation %s_%s_%s took %f seconds."
            % (wild, resNum, mutation, finishtime - starttime)
        )
        result = rosetta_binder.read_rosetta_ddgout(
            'mtfile.ddg', wild, mutation, resNum)
        # print(cartddg_cmd)
        os.chdir("../../")
        return result

    @staticmethod
    def run_row1(varlist: list):
        wild, mutation, resNum, jobID, relaxedpdb, exe, rosettadb = varlist
        path_job_id = common.ROSETTA_JOBS_DIR + jobID
        # print(f"path_job_id: {path_job_id}")
        distutils.dir_util.mkpath(path_job_id)
        os.chdir(path_job_id)
        # print(os.getcwd())
        os.system("cp ../../" + common.ROSETTA_RELAX_DIR + relaxedpdb + " ./")
        with open("mtfile", "w+") as mtfile:
            mtfile.write("total 1\n")
            mtfile.write("1\n")
            mtfile.write(wild + " " + str(resNum) + " " + mutation + "\n")
            mtfile.close()

        argument_list = [
            exe,
            "-database",
            rosettadb,
            "-use_input_sc",
            "-s",
            relaxedpdb,
            "-ddg::local_opt_only true",
            "-ddg::opt_radius 0.1",
            "-ddg::weight_file soft_rep_design",
            "-ddg::iterations 1",
            "-ddg::min_cst false",
            "-ddg::mean false",
            "-ddg:min true",
            "-ddg::sc_min_only false",
            "-ddg::mut_file mtfile",
            "1>/dev/null",
        ]
        # print(argument_list)
        ddg_cmd = " ".join(argument_list)
        # print(ddg_cmd)
        starttime = time.time()
        os.system(ddg_cmd)
        finishtime = time.time()
        print(
            "[DEBUG]: Rosetta mutation %s_%s_%s took %f seconds."
            % (wild, resNum, mutation, finishtime - starttime)
        )
        # result = rosetta_binder.read_rosetta_ddgout('mtfile.ddg', wild, mutation, resNum)
        result = rosetta_binder.read_ddg_monomer_out(
            'ddg_predictions.out', wild, mutation, resNum)
        os.chdir("../../")
        # print(f"After: {os.getcwd()}")
        return result


if __name__ == "__main__":
    print("run")
    # pdbname = '1PGA.pdb'
    # chain = 'A'
    # threads = 24
    # relax_num = 200
    # prot = Rosetta(pdbname, relax_num, 'numThreads')
    # score, std = prot.read_rosetta_ddgout('rosetta_jobs/0_1_/mtfile.ddg')
    # print(score, std)
    # relax_pdb_name = prot.relax(pdbname, threads)
    # print("Using: " + relax_pdb_name)
    # os.mkdir('rosetta_jobs')
    # prot.runOneJob(relax_pdb_name, "M", chain, "Q", 1, 'rosetta_jobs/0_1/')
