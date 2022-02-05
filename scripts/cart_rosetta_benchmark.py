#!/usr/bin/env python

def hash_rosettaRes_PdbRes(pdb, chain):
    '''Most PDB file adopted the biological residue numbering, 
    Rosetta using a more numerical numbering begin from 1
    input: a opened pdbfile or a pdbfile name
    output: hashmap(bio_res, index)'''
    resNumList = list()
    if type(pdb) == str:
        pdb = open(pdb)
    
    for line in pdb:
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

def rosetta_cart_relax(pdb, threads, nstruct):
    
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
    
    relax_cmd = " ".join(
        [
            "mpirun -n",
            str(int(threads)),
            "relax.mpi.linuxgccrelease -s",
            pdb,
            "-use_input_sc",
            "-constrain_relax_to_start_coords",
            "-ignore_unrecognized_res",
            "-nstruct",
            str(int(nstruct)),
            "-relax:coord_constrain_sidechains",
            "-relax:cartesian",
            "-score:weights ref2015_cart",
            "-relax:min_type lbfgs_armijo_nonmonotone",
            "-relax:script cart2.script 1>/dev/null && sort -nk2 score.sc |head -n 1|awk '{print$22}'",
        ]
    )
    return relax_cmd

def rosetta_cart_ddg(pdb):
    argument_list = [
        'cartesian_ddg.mpi.linuxgccrelease',
        "-database",
        '/opt/rosetta_bin_linux_2021.16.61629_bundle/main/database',
        "-use_input_sc",
        "-s",
        pdb,
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
    return cartddg_cmd

def read_big_mtfile(mtfile):
    all_mutations = []
    with open(mtfile) as mtfile:
        for line in mtfile:
            try:
                wt, pos, mut = line.replace("\n", "").split()
                all_mutations.append([wt, pos, mut])
            except ValueError:
                continue
        mtfile.close()
    return all_mutations

# 

def run_one_job(varlist):
    '''absolute path to pdb is required'''
    pdb, mutation, inv_res_dict = varlist
    wt, pos, mut = mutation
    
    bio_num = inv_res_dict[int(pos)]
    sub_dir = "_".join([wt, str(bio_num), mut])
    os.mkdir(sub_dir)
    # print("cp %s %s"%(pdb, sub_dir))
    os.popen("cp %s %s"%(pdb, sub_dir))
    os.chdir(sub_dir)
    
    with open("mtfile", 'w+') as mtfile:
        mtfile.write("total 1\n1\n%s %s %s\n"%(wt, str(pos), mut))
        mtfile.close()
        
    
    cart_ddg_cmd = rosetta_cart_ddg(pdb.split("/")[-1])
    # print(cart_ddg_cmd)
    starttime = time.time()
    os.system(cart_ddg_cmd)
    finishtime = time.time()
    ddg = read_rosetta_ddgout('mtfile.ddg')
            
    os.chdir("..")
    return [wt, bio_num, mut, ddg]
    
def run_all_jobs(pdb, mtfile, inv_res_dict, threads):
    all_results = []
    job_list = list()
    all_mutations = read_big_mtfile(mtfile)
    for mutation in all_mutations:
        job_list.append([pdb, mutation, inv_res_dict])

    # scan_start = time.time()
    all_results = Parallel(n_jobs=int(threads))(delayed(run_one_job)(var) for var in job_list)
    return all_results
    
def read_rosetta_ddgout(rosettaddgfilename):
    ddg_array = [[],[]]
    with open(rosettaddgfilename) as rosettaddg:
        for line in rosettaddg:
            # print(line.split(":")[2])
            if line.split(":")[2].strip() == "WT":
                dg_ref = float(line.split(":")[3][1:10])
                ddg_array[0].append(dg_ref)
            else:
                ddg = float(line.split(":")[3][1:10])
                ddg_array[1].append(ddg)
        rosettaddg.close()
    return np.array(ddg_array[1]).mean() - np.array(ddg_array[0]).mean()

def write_scores(all_results):
    with open("cart_scores.txt", 'w+') as scores:
        for x in all_results:
            scores.write("\t".join(x)+"\n")
        scores.close()

def main():
    pdb, chain, threads, nstruct = sys.argv[1:]
    res_dict = hash_rosettaRes_PdbRes(pdb, chain)
    inv_res_dict = {v: k for k, v in res_dict.items()} #get bio_numbering from rosetta_pos
    relax_cmd = rosetta_cart_relax(pdb, threads, nstruct)
    relaxed_pdb = os.popen(relax_cmd).read().strip() + ".pdb"
    # print(relaxed_pdb)
    # pdb = os.path.join(os.getcwd(),relaxed_pdb)
    all_results = run_all_jobs(relaxed_pdb, 'mtfile', inv_res_dict, threads)
    # print(all_results)
    write_scores(all_results)
    
if __name__ == '__main__':
    import sys
    import numpy as np
    import pandas as pd
    import time
    import os
    from joblib import Parallel, delayed
    # print()
    main()
