#!/usr/bin/env python

import pandas as pd
import os
import utlis.foldx as foldx
import utlis.rosetta as rosetta
from utlis import abacus
from utlis import mdrelax
import utlis.io as io
from joblib import Parallel, delayed


from multiprocessing import Pool, cpu_count


class GRAPE:
    def __init__(self):
        self.repaired_pdbfile: str
        self.relaxed_prot: str
        # self.repaired_pdbfile: str
        pass
    def run_foldx(self, pdb, threads, chain, numOfRuns):
        prot_foldx = foldx.FoldX(pdb, '', threads)
        self.repaired_pdbfile = prot_foldx.repairPDB()
        prot = io.Protein(self.repaired_pdbfile, chain)
        seq, resNumList = io.Protein.pdb2seq(prot)
        # parallel:
        # parall = ParallelSim(threads)
        try:
            os.mkdir('foldx_jobs')
        except FileExistsError:
            pass
        all_results = []
        job_list = []
        for i, res in enumerate(seq):
            resNum = resNumList[i]
            wild = res
            for j, aa in enumerate('QWERTYIPASDFGHKLCVNM'):
                jobID = "foldx_jobs/" + "_".join([wild, str(resNum), aa])
                job_list.append([self.repaired_pdbfile, wild, chain, aa, resNum, jobID, numOfRuns])

                # parall.add(prot_foldx.runOneJob, [self.repaired_pdbfile, wild, chain, aa, resNum, jobID])
        # parall.run()
        # all_results = parall.get_results()
        # print(job_list)
        Parallel(n_jobs=threads)(delayed(prot_foldx.runOneJob)(var) for var in job_list)
        
        return all_results

    def run_rosetta(self, pdb, threads, chain, relax_num):
        # relax_num = 200
        prot_rosetta = rosetta.Rosetta(pdb, relax_num, threads)
        relaxed_prot = prot_rosetta.relax()

        prot = io.Protein(pdb, chain)
        seq, resNumList = io.Protein.pdb2seq(prot)

        try:
            os.mkdir('rosetta_jobs')
        except FileExistsError:
            pass
        # all_results = []
        job_list = []
        for i, res in enumerate(seq):
            resNum = resNumList[i]
            wild = res
            for j, aa in enumerate('QWERTYIPASDFGHKLCVNM'):
                jobID = "rosetta_jobs/" + "_".join([wild, str(resNum), aa])
                job_list.append([wild, aa, str(i + 1), jobID])


        Parallel(n_jobs=threads)(delayed(prot_rosetta.runOneJob)(var) for var in job_list)

        return prot_rosetta

    def Analysis_foldx(self, pdb, chain, foldx1):
        self.repaired_pdbfile = pdb.replace(".pdb", '_Repair.pdb')
        try:
            os.mkdir("foldx_results")
        except FileExistsError:
            pass
        prot = io.Protein(pdb, chain)
        seq, resNumList = io.Protein.pdb2seq(prot)

        all_results = []
        for i, res in enumerate(seq):
            resNum = resNumList[i]
            wild = res
            for j, aa in enumerate('QWERTYIPASDFGHKLCVNM'):
                # jobID = "foldx_jobs/" + str(i) + "_" + str(j) + "/"
                jobID = "foldx_jobs/" + "_".join([wild, str(resNum), aa])
                all_results.append(foldx1.calScore(wild, resNum, aa, self.repaired_pdbfile, jobID))

        with open("foldx_results/All_FoldX.score", 'w+') as foldxout:
            foldxout.write("#Score file formated by GRAPE from FoldX.\n#mutation\tscore\tstd\n")
            for line in all_results:
                foldxout.write('\t'.join([line[0], str(line[1]), str(line[2])]) + "\n")
            foldxout.close()

        return all_results

    def Analysis_rosetta(self, pdb, chain, prot_rosetta):
        # self.repaired_pdbfile = pdb.replace(".pdb", '_Repair.pdb')
        # prot_rosetta = rosetta.Rosetta(pdb, relax_num, threads)
        try:
            os.mkdir("rosetta_results")
        except FileExistsError:
            pass
        prot = io.Protein(pdb, chain)
        seq, resNumList = io.Protein.pdb2seq(prot)

        all_results = []
        for i, res in enumerate(seq):
            resNum = resNumList[i]
            wild = res
            for j, aa in enumerate('QWERTYIPASDFGHKLCVNM'):
                # jobID = "foldx_jobs/" + str(i) + "_" + str(j) + "/"
                # "_".join([wild, str(resNum), mutation])
                rosettaddgfile = "rosetta_jobs/" + "_".join([wild, str(resNum), aa]) + "/mtfile.ddg"
                all_results.append(["_".join([wild, str(resNum), aa])] + prot_rosetta.read_rosetta_ddgout(rosettaddgfile))

        with open("rosetta_results/All_rosetta.score", 'w+') as rosettaout:
            rosettaout.write("#Score file formated by GRAPE from Rosetta.\n#mutation\tscore\tstd\n")
            for line in all_results:
                rosettaout.write('\t'.join([line[0], str(line[1]), str(line[2])]) + "\n")
            rosettaout.close()

        return all_results

    def analysisGrapeScore(self, scoreFile, cutoff, result_dir):
        result_dict = {'mutation':[],'energy':[],'SD':[],'position':[]}
        with open(scoreFile) as scorefile:
            for line in scorefile:
                if line[0] != "#":
                    lst = line.strip().split("\t")
                    result_dict['mutation'].append(lst[0].replace("_", ""))
                    result_dict['energy'].append(float(lst[1]))
                    result_dict['SD'].append(float(lst[2]))
                    result_dict['position'].append(int(lst[0].split("_")[1]))
            scorefile.close()
        #print(result_dict)
        CompleteList_df = pd.DataFrame(result_dict)
        CompleteList_SortedByEnergy_df = CompleteList_df.sort_values('energy').reset_index(drop=True)

        def BetsPerPosition(df):
            position_list = []
            length = df.shape[0]
            for i in range(length):
                if df['position'][i] in position_list:
                    df = df.drop(index=i)
                else:
                    position_list.append(df['position'][i])
            return df.reset_index(drop=True)

        def BelowCutOff(df,cutoff):
            #position_list = []
            length = df.shape[0]
            for i in range(length):
                if float(df['energy'][i]) > float(cutoff):
                    df = df.drop(index=i)
                else:
                    continue
            return df.reset_index(drop=True)

        BestPerPosition_SortedByEnergy_df = BetsPerPosition(CompleteList_SortedByEnergy_df)
        BestPerPosition_df = BetsPerPosition(CompleteList_SortedByEnergy_df)
        BelowCutOff_df = BelowCutOff(CompleteList_df, cutoff)
        BelowCutOff_SortedByEnergy_df = BelowCutOff(CompleteList_SortedByEnergy_df, cutoff)
        BestPerPositionBelowCutOff_SortedByEnergy_df = BelowCutOff(BestPerPosition_SortedByEnergy_df, cutoff)
        BestPerPositionBelowCutOff_df = BelowCutOff(BestPerPosition_df, cutoff)

        def out_tab_file(df, name, result_dir):
            filename = result_dir + "/MutationsEnergies_" + name[:-3] + ".tab"
            with open(filename,"w+") as of:
                of.write(df.to_csv(columns=['mutation', 'energy', 'SD'], sep='\t', index=False))
                of.close()

        out_tab_file(CompleteList_df, "CompleteList_df", result_dir)
        out_tab_file(CompleteList_SortedByEnergy_df, "CompleteList_SortedByEnergy_df", result_dir)
        out_tab_file(BestPerPosition_SortedByEnergy_df, "BestPerPosition_SortedByEnergy_df", result_dir)
        out_tab_file(BestPerPosition_df, "BestPerPosition_df", result_dir)
        out_tab_file(BelowCutOff_df, "BelowCutOff_df", result_dir)
        out_tab_file(BelowCutOff_SortedByEnergy_df, "BelowCutOff_SortedByEnergy_df", result_dir)
        out_tab_file(BestPerPositionBelowCutOff_SortedByEnergy_df, "BestPerPositionBelowCutOff_SortedByEnergy_df", result_dir)
        out_tab_file(BestPerPositionBelowCutOff_df, "BestPerPositionBelowCutOff_df", result_dir)

def selectpdb4md(pdb, platform, softlist):
    try:
        os.mkdir("selectpdb")
    except FileExistsError:
        os.system("rm -rf selectpdb")
    selected_dict = {"mutation": [], "score": [], "sd": [], "soft": []}
    for soft in softlist:
        with open("%s_results/MutationsEnergies_BelowCutOff.tab") as scorefile:
            for line in scorefile:
                linelist = line.strip().split()
                if linelist[0] != "mutation":
                    selected_dict["mutation"].append(linelist[0])
                    selected_dict["score"].append(linelist[1])
                    selected_dict["score"].append(linelist[2])
                    selected_dict["soft"].append(soft)
            scorefile.close()
    selected_df = pd.DataFrame(selected_dict)
    selected_df.to_csv("Selected_Mutation.csv")
    for mutation in selected_dict['mutation']:
        mut_pdb = pdb.replace(".pdb", "_Repair_1_0.pdb")
        os.system("cp foldx_jobs/%s/%s selectpdb/%s.pdb" %(mutation, mut_pdb, mutation))
        os.chdir("selectpdb")
        mutant = mutation + ".pdb"
        mdrelax.main(mutant, mutation + "_afterMD.pdb", platform)
        os.system("rm *dcd")
        os.chdir("..")


def main1():
    args = io.Parser().get_args()
    #print(args)

    pdb = args.pdb
    chain = args.chain
    threads = int(args.threads)
    numOfRuns = str(args.numofruns)
    # ratio = args.ratio
    relax_num = args.relax_number
    foldx_cutoff = -float(args.foldx_cutoff)
    rosetta_cutoff = -float(args.rosetta_cutoff)
    softlist = args.softlist.lower().split(",")
    preset = args.preset.lower()
    md = args.molecular_dynamics
    platform = args.platform


    exe_dict = {'foldx': '', 'relax': '', 'cartddg': '', 'pmut': '', 'abacus': ''}

    foldx_exe = os.popen("which foldx").read().replace("\n", "")
    exe_dict['foldx'] = foldx_exe
    pmut_scan_parallel_exe = os.popen("which pmut_scan_parallel.mpi.linuxgccrelease").read().replace("\n", "")
    rosettadb = "/".join(pmut_scan_parallel_exe.split("/")[:-3]) + "/database/"
    exe_dict['pmut'] = pmut_scan_parallel_exe
    for release in ['','.static','.mpi','.default']:
        cartesian_ddg_exe = os.popen("which cartesian_ddg%s.linuxgccrelease" %(release)).read().replace("\n", "")
        if cartesian_ddg_exe != "":
            exe_dict['cartddg'] = cartesian_ddg_exe
    relax_exe = os.popen("which relax.mpi.linuxgccrelease").read().replace("\n", "")
    exe_dict['relax'] = relax_exe
    abacus_prep = os.popen("which ABACUS_prepare").read().replace("\n", "")

    exe_dict['abacus'] = abacus_prep

    for soft in softlist:
        if soft == 'rosetta':
            if exe_dict['relax'] == '':
                print("[Error:] Cannot find Rosetta: relax.mpi.linuxgccrelease!")
                exit()
            if preset == 'slow':
                if exe_dict['cartddg'] == '':
                    print("[Error:] Cannot find Rosetta: any cartesian_ddg.linuxgccrelease (mpi nor default nor static)!")
                    exit()
            if preset == 'fast':
                if exe_dict['pmut'] == '':
                    print("[Error:] Cannot find Rosetta: pmut_scan_parallel.mpi.linuxgccrelease!")
                    exit()
        else:
            if exe_dict[soft] == '':
                print("[Error:] Cannot find %s!" %(soft))
                exit()



    mode = args.mode

    grape = GRAPE()
    foldx1 = foldx.FoldX(pdb, foldx_exe, threads)
    rosetta1 = rosetta.Rosetta(pdb, relax_num, threads, cartesian_ddg_exe, rosettadb)

    if mode == "run":
        #FoldX
        if "foldx" in softlist:
            grape.run_foldx(pdb, threads, chain, numOfRuns)
            grape.Analysis_foldx(pdb, chain, foldx1)
            grape.analysisGrapeScore('foldx_results/All_FoldX.score', foldx_cutoff, "foldx_results/")
        if preset == "slow":
            if "rosetta" in softlist:
                prot_rosetta = grape.run_rosetta(pdb, threads, chain, relax_num)
                grape.Analysis_rosetta(pdb, chain, prot_rosetta)
                grape.analysisGrapeScore('rosetta_results/All_rosetta.score', rosetta_cutoff, "rosetta_results/")
        if preset == "fast":
            if "rosetta" in softlist:
                relaxed_pdb = rosetta1.relax()

                try:
                    os.mkdir("rosetta_jobs")
                    os.chdir("rosetta_jobs")
                except FileExistsError:
                    os.chdir("rosetta_jobs")
                os.system("cp rosetta_relax/%s ./" %(relaxed_pdb))
                rosetta1.pmut_scan(relaxed_pdb)
                os.chdir("..")

                try:
                    os.mkdir("rosetta_results")
                    os.chdir("rosetta_results")
                except FileExistsError:
                    os.chdir("rosetta_results")
                    os.system("rm *")
                rosetta1.pmut_scan_analysis("../rosetta_jobs/pmut.out")
                os.chdir("..")

                grape.analysisGrapeScore('rosetta_results/All_rosetta.score', rosetta_cutoff, "rosetta_results/")

                # prot_rosetta = grape.run_rosetta(pdb, threads, chain, relax_num)
                # grape.Analysis_rosetta(pdb, chain, prot_rosetta)
                # grape.analysisGrapeScore('rosetta_results/All_rosetta.score', rosetta_cutoff, "rosetta_results/")
        if "abacus" in softlist:
            abacus.run_abacus(pdb)
            abacus.parse_abacus_out()
            grape.analysisGrapeScore('ABACUS_results/All_ABACUS.score', rosetta_cutoff, "ABACUS_results/")
    if mode == "analysis":
        #FoldX
        if "foldx" in softlist:
            # pdb = pdb.replace(".pdb", "_Repair.pdb")
            grape.Analysis_foldx(pdb, chain, foldx1)
            grape.analysisGrapeScore('foldx_results/All_FoldX.score', foldx_cutoff, "foldx_results/")
        if preset == "slow":
            if "rosetta" in softlist:
                prot_rosetta = rosetta.Rosetta(pdb, relax_num, threads)
                grape.Analysis_rosetta(pdb, chain, prot_rosetta)
                grape.analysisGrapeScore('rosetta_results/All_rosetta.score', rosetta_cutoff, "rosetta_results/")
        if preset == "fast":
            if "rosetta" in softlist:
                try:
                    os.mkdir("rosetta_results")
                    os.chdir("rosetta_results")
                except FileExistsError:
                    os.chdir("rosetta_results")
                    os.system("rm *")
                rosetta1.pmut_scan_analysis("../rosetta_jobs/pmut.out")
                os.chdir("..")
                grape.analysisGrapeScore('rosetta_results/All_rosetta.score', rosetta_cutoff, "rosetta_results/")
        if "abacus" in softlist:
            # abacus.run_abacus(pdb)
            abacus.parse_abacus_out()
            grape.analysisGrapeScore('abacus_results/All_ABACUS.score', rosetta_cutoff, "abacus_results/")


    print('\nFinish calculation of single point mutation ddG.\n')

    if md == True:
        selectpdb4md(pdb, platform, softlist)
    else:
        print("No MDs!")




if __name__ == '__main__':
    main1()
