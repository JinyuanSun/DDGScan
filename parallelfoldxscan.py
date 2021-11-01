#!/usr/bin/env python

import pandas as pd
import os
import utlis.foldx as foldx
import utlis.io as io

from multiprocessing import Pool, cpu_count


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

class GRAPE:
    def __init__(self):
        
        pass
    def run_foldx(self, pdb, threads, chain):
        prot_foldx = foldx.FoldX(pdb, '', threads, foldx_cutoff)
        pdbname = prot_foldx.repairPDB()
        prot = io.Protein(pdbname, chain)
        seq, resNumList = io.Protein.pdb2seq(prot)
        # parallel:
        parall = ParallelSim(threads)
        os.mkdir('foldx_jobs')
        all_results = []

        for i, res in enumerate(seq):
            resNum = resNumList[i]
            wild = res
            for j, aa in enumerate('QWERTYIPASDFGHKLCVNM'):
                jobID = "foldx_jobs/" + str(i) + "_" + str(j)
                parall.add(prot_foldx.runOneJob, [pdbname, wild, chain, aa, resNum, jobID])
        parall.run()
        all_results = parall.get_results()
        
        return all_results

    def Analysis_foldx(self, pdbname, chain, foldx1):
        try:
            os.mkdir("foldx_results")
        except FileExistsError:
            pass
        prot = io.Protein(pdbname, chain)
        seq, resNumList = io.Protein.pdb2seq(prot)

        all_results = []
        for i, res in enumerate(seq):
            resNum = resNumList[i]
            wild = res
            for j, aa in enumerate('QWERTYIPASDFGHKLCVNM'):
                jobID = "foldx_jobs/" + str(i) + "_" + str(j) + "/"
                all_results.append(foldx1.calScore(wild, resNum, aa, pdbname, jobID))
        
        with open("foldx_results/All_FoldX.score", 'w+') as foldxout:
            foldxout.write("#Score file formated by GRAPE.\n#mutation\tscore\tstd\n")
            for line in all_results:
                foldxout.write('\t'.join([line[0], str(line[1]), str(line[2])]) + "\n")
            foldxout.close()

        return all_results
    
    def analysisGrapeScore(self, scoreFile, cutoff):
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

        def out_tab_file(df, name):
            filename = "foldx_results/MutationsEnergies_" + name[:-3] + ".tab"
            with open(filename,"w+") as of:
                of.write(df.to_csv(columns=['mutation', 'energy', 'SD'], sep='\t', index=False))
                of.close()
                
        out_tab_file(CompleteList_df, "CompleteList_df")
        out_tab_file(CompleteList_SortedByEnergy_df, "CompleteList_SortedByEnergy_df")
        out_tab_file(BestPerPosition_SortedByEnergy_df, "BestPerPosition_SortedByEnergy_df")
        out_tab_file(BestPerPosition_df, "BestPerPosition_df")
        out_tab_file(BelowCutOff_df, "BelowCutOff_df")
        out_tab_file(BelowCutOff_SortedByEnergy_df, "BelowCutOff_SortedByEnergy_df")
        out_tab_file(BestPerPositionBelowCutOff_SortedByEnergy_df, "BestPerPositionBelowCutOff_SortedByEnergy_df")
        out_tab_file(BestPerPositionBelowCutOff_df, "BestPerPositionBelowCutOff_df")

if __name__ == '__main__':
    args = io.Parser.get_args()
    #print(args)
    
    pdb = args.pdb
    chain = args.chain
    threads = int(args.threads)
    # ratio = args.ratio
    foldx_cutoff = -float(args.foldx_cutoff)

    mode = args.mode

    grape = GRAPE()
    foldx1 = foldx.FoldX(pdb, '', threads, foldx_cutoff)


    if mode == "run":
        #FoldX
        grape.run_foldx(pdb, threads, chain)
        grape.Analysis_foldx(pdb, chain, foldx1)
        grape.analysisGrapeScore('All_FoldX.score', foldx_cutoff)
    if mode == "analysis":
        #FoldX
        pdb = pdb.replace(".pdb", "_Repair.pdb")
        grape.Analysis_foldx(pdb, chain, foldx1)
        grape.analysisGrapeScore('foldx_results/All_FoldX.score', foldx_cutoff)


    print('Done')




