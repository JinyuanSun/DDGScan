#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/3/28 11:25 PM
# @File    : list_distribute.py
# @annotation    :
# TODO: select mutation by relative properties

import argparse
import os

import pandas as pd
from joblib import Parallel, delayed

from utlis.foldx import foldx_binder
from utlis.rosetta import rosetta_binder

# from utlis import modeller_loop

class_type_dict = {
    '_small': 'GAVSTC',
    '_large': 'FYWKRHQE',
    '_neg': 'DE',
    '_pos': 'RK',
    '_polar': 'YTSHKREDQN',
    '_non_charged_polar': 'YTSNQH',
    '_hydrophobic': 'FILVAGMW',
    '_cys': "C",
    '_pro': 'P',
    '_scan': 'ARNDCQEGHILKMFPSTWYV'
}


def read_list(mutation_list_file):
    """
    param: mutation_list_file: a space-seperated text file
           wildtype chain position mutation
           A A 26 P
           A A 26 ILV
           A A 26 _polar
    return mutation_list
    """
    mutation_list = []
    with open(mutation_list_file, 'r') as mutation_list_file:
        for line in mutation_list_file:
            wildtype, chain, position, mutations = line.replace("\n", "").split(" ")
            if mutations[0] == "_":
                mutations = class_type_dict[mutations]
            for _, aa in enumerate(list(mutations)):
                if aa != wildtype:
                    mutation_list.append("_".join([wildtype,
                                                   chain,
                                                   position,
                                                   aa]))
        mutation_list_file.close()
        mutation_list = list(set(mutation_list))
        return mutation_list


class FoldX:
    def __init__(self):
        pass

    @staticmethod
    def mk_job_list(pdb_file, numOfRuns, mutation_list):
        job_list = []
        for mutation in mutation_list:
            wild, chain, position, mutation = mutation.split("_")
            job_id = "_".join([wild, position, mutation])
            var_list = [pdb_file, wild, chain, mutation, position, job_id, numOfRuns]
            # pdb_file, wild, chain, mutation, position, job_id, numOfRuns = varlist
            job_list.append(var_list)
        return job_list

    @staticmethod
    def dump_score_file(results, pdb):
        pdb_id = pdb.replace(".pdb", "")
        with open(pdb_id + "_FoldX.score", 'w+') as outfile:
            outfile.write('\t'.join(['mutation', 'mean', 'min', 'std']) + '\n')
            for index, result in enumerate(results):
                outfile.write("\t".join(result) + '\n')
            outfile.close()


class Rosetta:
    def __init__(self):
        pass

    @staticmethod
    def hash_rosettaRes_PdbRes(pdb, chain):
        '''Most PDB file adopted the biological residue numbering,
        Rosetta using a more numerical numbering begin from 1
        input: a opened pdbfile or a pdbfile name
        output: hashmap(bio_res, index)'''
        resNumList = list()
        if type(pdb) == str:
            pdb = open(pdb)

        for line in pdb:
            if "ATOM" == line[0:6].replace(" ", ""):
                if chain == line[21].replace(" ", ""):
                    if line[12:16].replace(" ", "") == "CA":
                        if line[16] == "B":
                            # print(line)
                            continue
                        else:
                            resNumList.append(int(line[22:26].replace(" ", "")))
        pdb.close()

        res_dict = {}
        for i, res in enumerate(resNumList):
            res_dict[res] = i + 1

        return res_dict

    @staticmethod
    def mk_job_list(pdb, relaxedpdb, mutation_list):

        def get_exe_db():
            relax_exe = os.popen("which relax.mpi.linuxgccrelease").read().replace("\n", "")
            rosettadb = os.popen("echo $ROSETTADB").read().replace("\n", "")
            if not rosettadb:
                rosettadb = "/".join(relax_exe.split("/")[:-4]) + "/database/"
            for release in ["", ".static", ".mpi", ".default"]:
                cartesian_ddg_exe = (
                    os.popen("which cartesian_ddg%s.linuxgccrelease" % (release))
                        .read()
                        .replace("\n", "")
                )
                if cartesian_ddg_exe != "":
                    # exe_dict["cartddg"] = cartesian_ddg_exe
                    return cartesian_ddg_exe, rosettadb

        exe, rosettadb = get_exe_db()

        job_list = []
        for mutation in mutation_list:
            wild, chain, position, mutation = mutation.split("_")
            res_dict = Rosetta.hash_rosettaRes_PdbRes(pdb, chain)
            resNum = res_dict[int(position)]
            job_id = "_".join([wild, position, mutation])
            var_list = [wild, mutation, resNum, job_id, relaxedpdb, exe, rosettadb]
            # wild, mutation, resNum, jobID, relaxedpdb, exe, rosettadb
            job_list.append(var_list)
        return job_list

    @staticmethod
    def dump_score_file(results, pdb):
        pdb_id = pdb.replace(".pdb", "")
        with open(pdb_id + "_Rosetta.score", 'w+') as outfile:
            outfile.write('\t'.join(['mutation', 'mean', 'min', 'std']) + '\n')
            for index, result in enumerate(results):
                outfile.write("\t".join(result) + '\n')
            outfile.close()


def get_args():
    parser = argparse.ArgumentParser(
        description="Run FoldX, Rosetta and ABACUS for in silico deep mutation scan."
    )
    parser.add_argument("pdb", help="Input PDB")
    parser.add_argument('mutation_list_file', help='Mutation list file, see README for details')
    parser.add_argument(
        "-msaddg",
        "--output_of_MSAddg",
        help="The format of MSAddg *.scan.txt, and there may be mismatch between your pdb and sequence",
        action='store_true',
    )
    parser.add_argument(
        "-fill",
        "--fill_break_in_pdb",
        help="Use modeller to fill missing residues in your pdb file. Use this option with caution!",
        action="store_true",
    )
    parser.add_argument(
        "-fix_mm",
        "--fix_mainchain_missing",
        help="fixing missing backbone bone using pdbfixer",
        action='store_true',
    )
    parser.add_argument(
        "-T",
        "--threads",
        help="Number of threads to run FoldX, Rosetta or ABACUS2",
        default=16,
        type=int,
    )
    parser.add_argument(
        "-nstruct",
        "--relax_number",
        help="Number of how many relaxed structure",
        default=50,
        type=int,
    )
    parser.add_argument(
        "-nruns",
        "--numofruns",
        help="Number of runs in FoldX BuildModel",
        default=5,
        type=int,
    )
    parser.add_argument(
        "-E",
        "--engine",
        nargs="+",
        choices=["foldx", "rosetta", "abacus2"],
    )
    parser.add_argument(
        "-repair",
        "--foldx_repair",
        help="Run Repair before ddG calculation",
        action="store_true",
    )
    parser.add_argument(
        "-MD",
        "--molecular_dynamics",
        help="Run 1ns molecular dynamics simulations for each mutation using openmm.",
        action="store_true",
    )
    parser.add_argument(
        "-P",
        "--platform",
        help="CUDA or CPU",
        type=str,
        choices=["CUDA", "CPU"],
        default="CUDA",
    )

    args = parser.parse_args()

    return args


def read_msaddg(msaddg_out, top=80, chain='A'):
    mutation_list = []
    df = pd.read_csv(msaddg_out, sep='\t', header=0, index_col=None)
    selected = df.sort_values('score', ascending=False).head(top)
    for mutation in selected['mutation'].values:
        wildtype, position, mutation = mutation
        mutation_list.append("_".join([wildtype, chain, position, mutation]))
    return mutation_list


def main(args=None):
    if not args:
        args = get_args()
    threads = args.threads
    pdb_file = args.pdb
    numOfRuns = args.numofruns
    repair = args.foldx_repair
    mutation_list_file = args.mutation_list_file
    output_of_MSAddg = args.output_of_MSAddg
    engines = args.engine
    relax_num = args.relax_number
    if output_of_MSAddg:
        mutation_list = output_of_MSAddg(mutation_list_file)
    else:
        mutation_list = read_list(mutation_list_file)
    if 'foldx' in engines:
        if repair:
            pdb_file = foldx_binder.repair_pdb(pdb_file)
        job_list = FoldX.mk_job_list(pdb_file, numOfRuns, mutation_list)
        results = Parallel(n_jobs=threads)(delayed(foldx_binder.run_one_job)(var) for var in job_list)
        FoldX.dump_score_file(results, args.pdb)
    if 'rosetta' in engines:
        relaxed_pdb = rosetta_binder.relax(args.pdb, threads, relax_num)
        job_list = Rosetta.mk_job_list(args.pdb, relaxed_pdb, mutation_list)
        results = Parallel(n_jobs=threads)(delayed(rosetta_binder.run_one_job)(var) for var in job_list)
        Rosetta.dump_score_file(results, args.pdb)


if __name__ == '__main__':
    main()
