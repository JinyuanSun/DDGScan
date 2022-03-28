#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/3/28 11:25 PM
# @File    : list_distribute.py
# @annotation    :
# TODO: select mutation by relative properties

import os
# from utlis import modeller_loop
import time

import pandas as pd
from joblib import Parallel, delayed
import argparse
from utlis.foldx import foldx_binder


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


def mk_job_list(pdb_file, numOfRuns, mutation_list):
    job_list = []
    for mutation in mutation_list:
        wild, chain, position, mutation = mutation.split("_")
        job_id = "_".join([wild, position, mutation])
        var_list = [pdb_file, wild, chain, mutation, position, job_id, numOfRuns]
        # pdb_file, wild, chain, mutation, position, job_id, numOfRuns = varlist
        job_list.append(var_list)
    return job_list


def dump_score_file(results):
    for index, result in enumerate(results):
        print("\t".join(result))


def get_args():
    parser = argparse.ArgumentParser(
        description="Run FoldX, Rosetta and ABACUS for in silico deep mutation scan."
    )
    parser.add_argument("pdb", help="Input PDB")
    parser.add_argument('mutation_list_file', help='Mutation list file, see README for details')
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





if __name__ == '__main__':
    args = get_args()
    threads = args.threads
    pdb_file = args.pdb
    numOfRuns = args.numofruns
    repair = args.foldx_repair
    mutation_list_file = args.mutation_list_file
    # threads = 8
    # numOfRuns = 5
    # pdb_file = '1pga.pdb'
    mutation_list = read_list(mutation_list_file)
    job_list = mk_job_list(pdb_file, numOfRuns, mutation_list)
    if repair:
        pdb_file = foldx_binder.repair_pdb(pdb_file)
    results = Parallel(n_jobs=threads)(delayed(foldx_binder.run_one_job)(var) for var in job_list)
    print(results)
    dump_score_file(results)
