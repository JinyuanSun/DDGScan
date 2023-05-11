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
from shutil import which

from utils.foldx import foldx_binder
from utils.rosetta import rosetta_binder
from utils.aa_index import *
from utils.common import ABACUS2_JOBS_DIR, ROSETTA_RELAX_DIR
import utils.abacus as abacus
from utils.abacus2_nn import *


# from utils import modeller_loop
from Bio.PDB import PDBParser, Select, PDBIO
from Bio.PDB.Polypeptide import one_to_three


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
    for aa, value in hydrophobic_index.items():
        if aa != wildtype:
            if mutation_type == '@less_hydrophobic':
                if value <= wt_hydropathy:
                    mutations += aa
            if mutation_type == '@more_hydrophobic':
                if value >= wt_hydropathy:
                    mutations += aa

    wt_sheet_tendenvy = sheet_tendency[wildtype]
    for aa, value in sheet_tendency.items():
        if aa != wildtype:
            if mutation_type == '@less_sheet_tendency':
                if value <= wt_sheet_tendenvy:
                    mutations += aa
            if mutation_type == '@more_sheet_tendency':
                if value >= wt_sheet_tendenvy:
                    mutations += aa

    wt_helix_tendency = helix_tendency[wildtype]
    for aa, value in helix_tendency.items():
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

def chain_resseq_to_pos_number(pdb):
    res_key_dict = {}
    parser = PDBParser(QUIET=True)
    model_0 = parser.get_structure("x", pdb)[0]
    pos_num = 0
    for chain in model_0.get_chains():
        chain_id = chain.id
        for residue in chain.get_residues():
            if residue.id[0] == " ":
                key = f"{chain_id}_{residue.id[1]}"
                pos_num += 1
                res_key_dict[key] = pos_num
    return res_key_dict

class Rosetta:
    def __init__(self):
        pass

    @staticmethod
    def hash_rosettaRes_PdbRes(pdb, chain):
        '''Most PDB file adopted the biological residue numbering,
        Rosetta using a more numerical numbering begin from 1
        input: a opened pdbfile or a pdbfile name
        output: hashmap(chain_resseq, index)'''
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
    def mk_job_list(pdb, relaxedpdb, mutation_list, fast=False):

        res_key_dict = chain_resseq_to_pos_number(pdb)

        def get_exe_db():
            relax_exe = which("relax.mpi.linuxgccrelease")
            rosettadb = os.popen("echo $ROSETTADB").read().replace("\n", "")
            if not rosettadb:
                rosettadb = "/".join(relax_exe.split("/")[:-4]) + "/database/"
            for release in ["", ".static", ".mpi", ".default"]:
                cartesian_ddg_exe = which(f"cartesian_ddg{release}.linuxgccrelease")
                if cartesian_ddg_exe != None:
                    return cartesian_ddg_exe, rosettadb

        exe, rosettadb = get_exe_db()
        if fast:
            for release in ["", ".static", ".mpi", ".default"]:
                exe = which(f"ddg_monomer{release}.linuxgccrelease")
                if exe != None:
                    break

        job_list = []
        for mutation in mutation_list:
            wild, chain, position, mutation = mutation.split("_")
            # res_dict = Rosetta.hash_rosettaRes_PdbRes(pdb, chain)
            resNum = res_key_dict[f"{chain}_{int(position)}"]
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
        choices=["foldx", "rosetta", "abacus2", "rosetta_fast"],
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

def mk_abacus2nn_joblist(pdb_file, mutation_list):
    job_list = []
    mutations = []
    for mutation in mutation_list:
        wild, chain, position, mutation = mutation.split("_")
        mutations.append(f"{wild}_{position}_{mutation}")
        var_list = [pdb_file, chain, position, one_to_three(mutation)]
        # pdb_file, wild, chain, mutation, position, job_id, numOfRuns = varlist
        job_list.append(var_list)
    return job_list, mutations

def dump_abacus_score_file(results, pdb):
    pdb_id = pdb.replace(".pdb", "")
    with open(pdb_id + "_ABACUS2.score", 'w+') as outfile:
        outfile.write('\t'.join(['mutation', 'total_score']) + '\n')
        for index, result in enumerate(results):
            outfile.write(f"{result[0]}\t{result[1]}\n")
        outfile.close()

def dump_abacus2nn_score_file(abacus2_results, pdb:str, mutations:list):
    pdb_id = pdb.replace(".pdb", "")
    df = pd.DataFrame(abacus2_results)
    df[5] = mutations
    df[['mutations', 'sai', 's1', 's2','pack', 'hb']] = df[[5,0,1,2,3,4]]
    df[['mutations', 'sai', 's1', 's2','pack', 'hb']].to_csv(os.path.join(ABACUS2_JOBS_DIR, 'abacus2_raw.csv'), sep=',', index=None)
    model_list = get_models()
    with torch.no_grad():
        all_preds = []
        for net in model_list:
            x = torch.tensor(np.array(abacus2_results)).float()
            pred0 = net(x)
            pred0 = pred0.ravel().numpy()
            all_preds.append(pred0)
    avg_pred = np.mean(all_preds, axis=0) + 0.6169168
    std_pred = np.std(all_preds, axis=0)
    min_pred = np.min(all_preds, axis=0) + 0.6169168
    with open(pdb_id + "_ABACUS2_NN.score", 'w+') as outfile:
            outfile.write('\t'.join(['mutation', 'mean', 'min', 'std']) + '\n')
            for mutation, avg_v, min_v, std in zip(mutations, avg_pred, min_pred, std_pred):
                outfile.write("{:}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(mutation, avg_v, min_v, std))
            outfile.close()


class ProSelect(Select):
    def accept_residue(self, residue):
        if residue.id[0] == ' ':
            return 1
        else:
            return 0

def clean_pdb(pdb):
    parser = PDBParser(QUIET=True)
    model_0 = parser.get_structure('x', pdb)[0]
    io=PDBIO()
    io.set_structure(model_0)
    oname = pdb.replace(".pdb", "_protein.pdb")
    io.save(oname, ProSelect())
    print("PDB cleaned.")
    return oname

def main(args):
    # args = get_args()
    threads = args.threads
    pdb_file = args.pdb
    numOfRuns = args.numofruns
    repair = args.foldx_repair
    mutation_list_file = args.mutation_list_file
    output_of_MSAddg = args.output_of_MSAddg
    engines = args.engine
    relax_num = args.relax_number
    relax = args.rosetta_relax
    pdb_file = clean_pdb(pdb_file)
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
        relaxed_pdb = rosetta_binder.relax(pdb_file, threads, relax_num)
        job_list = Rosetta.mk_job_list(pdb_file, relaxed_pdb, mutation_list)
        results = Parallel(n_jobs=threads)(delayed(rosetta_binder.run_one_job)(var) for var in job_list)
        Rosetta.dump_score_file(results, args.pdb)
    if 'rosetta_fast' in engines:
        if relax:
            pdb_file = rosetta_binder.fast_relax(pdb_file, threads, relax_num)
        else:
            distutils.dir_util.mkpath(ROSETTA_RELAX_DIR)
            os.system(f"cp  {pdb_file} {ROSETTA_RELAX_DIR}")
        job_list = Rosetta.mk_job_list(pdb_file, pdb_file, mutation_list, fast=True)
        results = Parallel(n_jobs=threads)(delayed(rosetta_binder.run_row1)(var) for var in job_list)
        Rosetta.dump_score_file(results, args.pdb)
    if 'abacus2' in engines:
        distutils.dir_util.mkpath(ABACUS2_JOBS_DIR)
        # abacus.runOneJob
        job_list = mk_abacus_joblist(pdb_file, mutation_list)
        abacus2_results = Parallel(n_jobs=threads)(delayed(abacus.runOneJob)(var) for var in job_list)
        dump_abacus_score_file(abacus2_results, args.pdb)
    if 'abacus2_nn' in engines:
        distutils.dir_util.mkpath(ABACUS2_JOBS_DIR)
        job_list, mutations = mk_abacus2nn_joblist(pdb_file, mutation_list)
        abacus2_results = Parallel(n_jobs=threads)(delayed(abacus.run_abacus2_cmd)(*var) for var in job_list)
        dump_abacus2nn_score_file(abacus2_results, args.pdb, mutations)



if __name__ == '__main__':
    args = get_args()
    main(args)
