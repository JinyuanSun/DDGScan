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
           A A 12 P
           A A 12 ILV
           A A 13 _polar
    return mutation_list
    """
    mutation_list = []
    with open(mutation_list_file, 'r') as mutation_list_file:
        for line in mutation_list_file:
            wildtype, chain, position, mutation = line.replace("\n", "").split(" ")
            if mutation[0] == "_":
                mutations = class_type_dict[mutation]
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


if __name__ == '__main__':
    threads = 8
    numOfRuns = 5
    mutation_list_file = 'test_file.txt'
    pdb_file = '../test/1pga.pdb'
    mutation_list = read_list(mutation_list_file)
    job_list = mk_job_list(pdb_file, numOfRuns, mutation_list)
    repair = False
    if repair:
        pdb_file = foldx_binder.repair_pdb(pdb_file)
    results = Parallel(n_jobs=threads)(delayed(foldx_binder.run_one_job)(var) for var in job_list)
    dump_score_file(results)
