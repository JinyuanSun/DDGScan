#!/usr/bin/env python

import pandas as pd
import os
import utlis.foldx as foldx
import utlis.rosetta as rosetta
import utlis.io as io
from joblib import Parallel, delayed

from multiprocessing import Pool, cpu_count


'''
class ParallelSim(object):
    def __init__(self, processes=cpu_count()):
        self.pool = Pool(processes=processes)
        self.total_processes = 0
        self.completed_processes = 0
        self.results = []

    def add(self, func, args):
        self.pool.apply_async(func=func, args=args, callback=self.complete)
        self.total_processes += 1

    def complete(self, result):
        self.results.extend(result)
        self.completed_processes += 1
        print('Progress: {:.2f}%'.format((self.completed_processes/self.total_processes)*100))

    def run(self):
        self.pool.close()
        self.pool.join()

    def get_results(self):
        return self.results
'''

class GRAPE:
    def __init__(self):
        self.repaired_pdbfile: str
        self.relaxed_prot: str
        # self.repaired_pdbfile: str
        pass
    def run_foldx(self, pdb, threads, chain, numOfRuns):
        prot_foldx = foldx.FoldX(pdb, '', threads, foldx_cutoff)
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

if __name__ == '__main__':
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
    softlist = args.softlist.split(",")

    mode = args.mode

    grape = GRAPE()
    foldx1 = foldx.FoldX(pdb, '', threads, foldx_cutoff)
    rosetta1 = rosetta.Rosetta(pdb, relax_num, threads)


    if mode == "run":
        #FoldX
        if "FoldX" in softlist:
            grape.run_foldx(pdb, threads, chain, numOfRuns)
            grape.Analysis_foldx(pdb, chain, foldx1)
            grape.analysisGrapeScore('foldx_results/All_FoldX.score', foldx_cutoff, "foldx_results/")
        if "Rosetta" in softlist:
            prot_rosetta = grape.run_rosetta(pdb, threads, chain, relax_num)
            grape.Analysis_rosetta(pdb, chain, prot_rosetta)
            grape.analysisGrapeScore('rosetta_results/All_rosetta', rosetta_cutoff, "rosetta_results/")
    if mode == "analysis":
        #FoldX
        if "FoldX" in softlist:
            pdb = pdb.replace(".pdb", "_Repair.pdb")
            grape.Analysis_foldx(pdb, chain, foldx1)
            grape.analysisGrapeScore('foldx_results/All_FoldX.score', foldx_cutoff)
        if "Rosetta" in softlist:
            prot_rosetta = rosetta.Rosetta(pdb, relax_num, threads)
            grape.Analysis_rosetta(pdb, chain, prot_rosetta)
            grape.analysisGrapeScore('rosetta_results/All_rosetta.score', rosetta_cutoff)


    print('Done')
