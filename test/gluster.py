#!/usr/bin/env python
# 22 Oct 2021
# Clustering mutations according to dTm, location and interactions for GRAPE

import os
import sys

import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.cluster import KMeans


class Protein:
    def __init__(self, pdbname, chain):
        self.pdbname = pdbname
        self.chain = chain
        self.seq, self.resNumList = self.pdb2seq()

    def _3_2_1(self, x):
        d = {
            "CYS": "C",
            "ASP": "D",
            "SER": "S",
            "GLN": "Q",
            "LYS": "K",
            "ILE": "I",
            "PRO": "P",
            "THR": "T",
            "PHE": "F",
            "ASN": "N",
            "GLY": "G",
            "HIS": "H",
            "LEU": "L",
            "ARG": "R",
            "TRP": "W",
            "ALA": "A",
            "VAL": "V",
            "GLU": "E",
            "TYR": "Y",
            "MET": "M",
        }
        y = d[x]
        return y

    def pdb2seq(self):
        # return a sequence and a resnumlist
        seq = ""
        resNumList = []
        with open(self.pdbname) as pdbfile:
            for line in pdbfile:
                if "ATOM" == line[0:6].replace(" ", ""):
                    if self.chain == line[21].replace(" ", ""):
                        if line[12:16].replace(" ", "") == "CA":
                            if line[16] == "B":
                                continue
                            else:
                                seq += self._3_2_1(line[17:20].replace(" ", ""))
                                # x = float(line[30:38].replace(" ", ""))
                                # y = float(line[38:46].replace(" ", ""))
                                # z = float(line[46:54].replace(" ", ""))
                                resNumList.append(int(line[22:26].replace(" ", "")))
        pdbfile.close()
        return seq, resNumList

    def readPdbCaCoord(self):
        CaCoordDict = {}
        with open(self.pdbname) as pdbfile:
            for line in pdbfile:
                if "ATOM" == line[0:6].replace(" ", ""):
                    if self.chain == line[21].replace(" ", ""):
                        if line[12:16].replace(" ", "") == "CA":
                            res = self._3_2_1(line[17:20].replace(" ", "")) + line[
                                22:26
                            ].replace(" ", "")
                            x = float(line[30:38].replace(" ", ""))
                            y = float(line[38:46].replace(" ", ""))
                            z = float(line[46:54].replace(" ", ""))
                            CaCoordDict[res] = np.array([x, y, z])
        pdbfile.close()
        return CaCoordDict


class Cluster:
    def __init__(self, clusterInput, pdbFile, chain, clusterNum, path2foldx, yasaraExe):
        self.clusterInput = clusterInput
        self.path2foldx = path2foldx
        self.yasaraExe = yasaraExe
        self.pdbFile = pdbFile
        self.chain = chain
        self.clusterNum = int(clusterNum)
        self.dataDic = {
            "mutation": [],
            "dTm": [],
            "Hbond": [],
            "hydrophobic": [],
            "entropy": [],
            "CA_x": [],
            "CA_y": [],
            "CA_z": [],
        }

    def readClusInput(self, clusterInput):
        # Args: clusterInput is a tab-separated file of mutation name and measured dTm
        # Out: a 2D list of [[wildtype residue, Reside Number in pdb, mutated residue, dtm]]
        positionList = []
        with open(clusterInput) as clusin:
            for line in clusin:
                mutation, dtm = line.strip().split(",")
                wild = mutation[0]
                mut = mutation[-1]
                resNum = int(mutation[1:-1])
                dtm = float(dtm)
                positionList.append([wild, resNum, mut, dtm])
            clusin.close()
        return positionList

    def dumpmcr(self):
        with open("cluster.mcr", "w+") as mcrfile:
            mcrfile.write("LOADPDB (MacroTarget)\n")
            mcrfile.write("Clean ALL\n")
            mcrfile.write(
                "ListHBoRes res (MutatedResidue), protein with distance < 5 from res (MutatedResidue)\n"
            )
            mcrfile.write(
                "ListIntRes res (MutatedResidue), protein with distance < 5 from res (MutatedResidue), Type=hydrophobic\n"
            )
            mcrfile.write("Exit\n")
            mcrfile.close()

    def inputCheck(self, positionList, pdbFile, chain):
        # Check wether the input mutation matches with pdbfile
        prot = Protein(pdbFile, chain)
        seq, resNumList = prot.pdb2seq()
        CaCoorDict = prot.readPdbCaCoord()
        reslist = []
        for i in range(len(seq)):
            reslist.append(seq[i] + str(resNumList[i]))

        for mutation in positionList:
            res = mutation[0] + str(mutation[1])
            if res in reslist:
                continue
            else:
                print("[ERROR]: " + res + "Not Found in PDB!")
                return False

        return CaCoorDict

    def buildModel(self, path2foldx, positionList):
        foldx = FoldX(self.pdbFile, path2foldx)
        pdbFile = foldx.repairPDB()
        for mutation in positionList:
            jobID = "selectpdb/" + mutation[0] + str(mutation[1]) + mutation[2]
            foldx.runOneJob(
                pdbFile, mutation[0], self.chain, mutation[2], str(mutation[1]), jobID
            )

    def extractYasaraOut(self, yasaraOut):
        # to deal with nasty yasara output
        outList = yasaraOut.split("\n")
        HbondNum = 0
        hydrophobicNum = 0
        for line in outList:
            if "listed." in line:
                hydrophobicNum = int(line.split(" ")[0])
            if "hydrogen bond" in line:
                HbondNum = int(line.split(" ")[0])
        return HbondNum, hydrophobicNum

    def GetHbondAndHydrophobic(self, mutation, pdbFile):
        # mutHbond - wildHbond, mutHydrophobic - wildHydrophobic
        # gla_af_clean_Repair.pdb -> gla_af_clean_Repair_N38LWOW.pdb
        pos = str(mutation[1])
        # wild = mutation[0] + pos
        # mut = mutation[0] + pos + mutation[2]
        mutfilename = "selectpdb/" + mutation[0] + pos + mutation[2] + ".pdb"
        # Yasara: /home/cuily/soft/yasara/yasara -txt cluster.mcr "MacroTarget = 'gla_af_clean.pdb'" "MutatedResidue = 19 "
        wildCmd = (
            '%s -txt cluster.mcr "MacroTarget = \'%s\'" "MutatedResidue = %s "'
            % (yasaraExe, pdbFile, pos)
        )

        wildOut = os.popen(wildCmd).read()
        wildHbond, wildHydrophobic = Cluster.extractYasaraOut(self, wildOut)

        mutCmd = '%s -txt cluster.mcr "MacroTarget = \'%s\'" "MutatedResidue = %s "' % (
            yasaraExe,
            mutfilename,
            pos,
        )
        print(wildCmd + "\n" + mutCmd)
        mutOut = os.popen(mutCmd).read()
        mutHbond, mutHydrophobic = Cluster.extractYasaraOut(self, mutOut)

        return mutHbond - wildHbond, mutHydrophobic - wildHydrophobic

    def entropy(self, mutationlist):
        if "P" in mutationlist:
            return 1
        elif "G" in mutationlist:
            return 1
        else:
            return 0

    def GenFeatures(self, positionList, CaCoorDict):

        for mutation in positionList:
            # position = str(mutation[1])
            res = mutation[0] + str(mutation[1])
            self.dataDic["mutation"].append(res + mutation[2])
            self.dataDic["dTm"].append(float(mutation[3]))
            print("[INFO]: Extract features of mutation %s" % ((res + mutation[2])))

            # get entropy
            # where is P or G, where is entropy
            # wild, resNum, mut, dtm
            self.dataDic["entropy"].append(self.entropy(mutation))

            # get coord
            coord = CaCoorDict[res]
            self.dataDic["CA_x"].append(coord[0])
            self.dataDic["CA_y"].append(coord[1])
            self.dataDic["CA_z"].append(coord[2])

            # get Hbond and hydrophobic contacts
            Hbond, Hydro = Cluster.GetHbondAndHydrophobic(
                self.yasaraExe, mutation, self.pdbFile
            )
            self.dataDic["Hbond"].append(Hbond)
            self.dataDic["hydrophobic"].append(Hydro)
        df = pd.DataFrame(self.dataDic)
        df.to_csv("features.csv", sep="\t", index=True, index_label="Index")
        return self.dataDic

    def clustering(self, featureFile, N_CLUSTERS):
        df = pd.read_csv(featureFile, sep="\t")
        labels = np.array(df["mutation"])
        data = np.array(
            df[["dTm", "Hbond", "hydrophobic", "entropy", "CA_x", "CA_y", "CA_z"]]
        )

        # data = preprocessing.normalize(data, norm='l2')
        # N_CLUSTERS = 3 #input

        kmeans = KMeans(init="k-means++", n_clusters=N_CLUSTERS, n_init=9)

        data = preprocessing.normalize(data, norm="l2")

        kmeans.fit(data)
        pred_classes = kmeans.predict(data)
        # Print out results
        with open("cluster.txt", "w+") as clusterout:
            for cluster in range(N_CLUSTERS):
                print("cluster: ", cluster)
                clusterout.write("cluster: " + str(cluster))
                print(labels[np.where(pred_classes == cluster)])
                clusterout.write(str(labels[np.where(pred_classes == cluster)]) + "\n")
            clusterout.close()

    def run(self):
        positionList = self.readClusInput(self.clusterInput)
        CaCoorDict = self.inputCheck(positionList, self.pdbFile, self.chain)
        self.buildModel(self.path2foldx, positionList)
        self.dumpmcr()
        dataDict = self.GenFeatures(positionList, CaCoorDict)
        self.clustering("features.csv", self.clusterNum)


class FoldX:
    def __init__(self, pdbName, path2foldx):
        self.pdbname = pdbName
        self.foldxExe = path2foldx
        self.result = []

    def repairPDB(self):
        cmd = self.foldxExe + " --command=RepairPDB --pdb=" + self.pdbname
        # print(cmd)
        print("==" * 20)
        print("Stage 1: Repair PDB file using FoldX")
        print("==" * 20)
        print("[INFO]: Running FoldX RepairPDB.")

        FoldX_out = os.popen(cmd).read()
        with open(".foldx_repair.log", "w+") as outfile:
            outfile.write(FoldX_out)
            outfile.close()

        print("[INFO]: FoldX RepairPDB Done!")
        os.mkdir("selectpdb")
        return self.pdbname.replace(".pdb", "_Repair.pdb")

    def runOneJob(self, pdbfile, wild, chain, mutation, resNum, jobID):
        # cwd: selectpdb/ cluster.txt upload.pdb Repaired.pdb
        #     /__selectpdb: mutation_foldx_jobs/ mutation.pdb(s)
        #                  /__: mutation(s)/
        os.mkdir(jobID)  # mkdir selectpdb/mutation
        os.chdir(jobID)  #
        with open("individual_list.txt", "w+") as indFile:
            indFile.write(wild + chain + str(resNum) + mutation + ";")
            indFile.close()
        cmd1 = "cp ../../" + pdbfile + " ./"
        mutPdbName = "".join([wild, str(resNum), mutation, ".pdb"])
        os.system(cmd1)
        cmd2 = (
            self.foldxExe
            + " --command=BuildModel --numberOfRuns=1 --mutant-file=individual_list.txt --pdb="
            + pdbfile
        )
        os.system(cmd2)
        os.system("cp " + pdbfile.replace(".pdb", "_1.pdb") + " ../" + mutPdbName)
        os.chdir("../../")
        return mutPdbName


if __name__ == "__main__":
    try:
        clusterInput, pdbFile, chain, clusterNum, foldxExe, yasaraExe = sys.argv[1:]
    except ValueError:
        try:
            clusterInput, pdbFile, chain = sys.argv[1:]
            clusterNum = 3
            foldxExe = "foldx"
            yasaraExe = "yasara"
        except ValueError:
            print(
                "Usage: python3 cluster.py <dtmData.txt> <protein.pdb> <Chain> <ClusterNum> <path/to/foldx> </path/to/yasara>"
            )
            print(
                "Warning: Your should at least provide <dtmData.txt> <protein.pdb> <Chain>"
            )
            print(
                "         If Yasara and FoldX both set in env, the default Cluster Number will be 3!"
            )
            exit()

    clus = Cluster(clusterInput, pdbFile, chain, clusterNum, foldxExe, yasaraExe)
    clus.run()
