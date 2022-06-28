#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/3/28 11:25 PM
# @File    : list_distribute.py
# @annotation    :

import argparse
import os

import distutils.dir_util

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from utils.foldx import foldx_binder
from utils.rosetta import rosetta_binder
from utils.aa_index import *
from utils.common import ABACUS2_JOBS_DIR
import utils.abacus as abacus

# from utils import modeller_loop


def convert_by_property_selection(wildtype, mutation_type):
    """
    @smaller: mutation to AA with smaller vdw
    @bigger: mutation to AA with bigger vdw
    @more_hydrophobic: mutation to AA more hydrophobic
    @less_hydrophobic: mutation to AA more hydrophilic
    @more_sheet_tendency: mutation to AA with higher sheet tendency
    @less_sheet_tendency: mutation to AA with higher sheet tendency
    @more_helix_tendency: mutation to AA with higher helix tendency
    @less_helix_tendency: mutation to AA with higher helix tendency
    @{random}: random is an integer in range 1 to 19 ,randomly select few mutations for you, good luck!
    param: one-letter token of amino acid
    return: list of amino acid
    """

    mutations = ''
    wt_volume = volume_index[wildtype]
    for aa, value in volume_index.items():
        if aa != wildtype:
            if mutation_type == '@smaller':
                if value <= wt_volume:
                    mutations += aa
            if mutation_type == '@bigger':
                if value >= wt_volume:
                    mutations += aa

    wt_hydropathy = hydrophobic_index[wildtype]
    for aa, value in hydrophobic_index.item():
        if aa != wildtype:
            if mutation_type == '@less_hydrophobic':
                if value <= wt_hydropathy:
                    mutations += aa
            if mutation_type == '@more_hydrophobic':
                if value >= wt_hydropathy:
                    mutations += aa

    wt_sheet_tendenvy = sheet_tendency[wildtype]
    for aa, value in sheet_tendency.item():
        if aa != wildtype:
            if mutation_type == '@less_sheet_tendency':
                if value <= wt_sheet_tendenvy:
                    mutations += aa
            if mutation_type == '@more_sheet_tendency':
                if value >= wt_sheet_tendenvy:
                    mutations += aa

    wt_helix_tendency = helix_tendency[wildtype]
    for aa, value in helix_tendency.item():
        if aa != wildtype:
            if mutation_type == '@less_helix_tendency':
                if value <= wt_helix_tendency:
                    mutations += aa
            if mutation_type == '@more_helix_tendency':
                if value >= wt_helix_tendency:
                    mutations += aa

    try:
        random_number = int(mutation_type.replace("@", ""))
        aa_list = ALPHABET.replace(wildtype, "")
        while random_number > 0:
            mutations += aa_list.pop(np.random.randint(len(aa_list)))
            random_number -= 1
    except ValueError:
        pass

    return mutations


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
            if mutations[0] == "@":
                mutations = convert_by_property_selection(wildtype, mutations)
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


def read_msaddg(msaddg_out, top=160, chain='A'):
    mutation_list = []
    df = pd.read_csv(msaddg_out, sep='\t', header=0, index_col=None)
    selected = df.sort_values('score', ascending=False).head(top)
    for mutation in selected['mutation'].values:
        wildtype, position, mutation = mutation.split('_')
        mutation_list.append("_".join([wildtype, chain, position, mutation]))
    return mutation_list

def mk_abacus_joblist(pdb_file, mutation_list):
    job_list = []
    for mutation in mutation_list:
        wild, chain, position, mutation = mutation.split("_")
        job_id = "_".join([wild, position, mutation])
        var_list = [pdb_file, wild, chain, mutation, position]
        # pdb_file, wild, chain, mutation, position, job_id, numOfRuns = varlist
        job_list.append(var_list)
    return job_list

def dump_abacus_score_file(results, pdb):
    pdb_id = pdb.replace(".pdb", "")
    with open(pdb_id + "_ABACUS2.score", 'w+') as outfile:
        outfile.write('\t'.join(['mutation', 'total_score']) + '\n')
        for index, result in enumerate(results):
            outfile.write(f"{result[0]}\t{result[1]}\n")
        outfile.close()


def main(args):
    threads = args.threads
    pdb_file = args.pdb
    numOfRuns = args.numofruns
    repair = args.foldx_repair
    mutation_list_file = args.mutation_list_file
    output_of_MSAddg = args.output_of_MSAddg
    engines = args.engine
    relax_num = args.relax_number
    if output_of_MSAddg:
        mutation_list = read_msaddg(mutation_list_file)
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

    if 'abacus2' in engines:
        distutils.dir_util.mkpath(ABACUS2_JOBS_DIR)
        # abacus.runOneJob
        job_list = mk_abacus_joblist(args.pdb, mutation_list)
        abacus2_results = Parallel(n_jobs=threads)(delayed(abacus.runOneJob)(var) for var in job_list)
        dump_abacus_score_file(abacus2_results, args.pdb)


if __name__ == '__main__':
    args = get_args()
    main(args)
