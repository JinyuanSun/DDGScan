#!/usr/bin/env python
import os
import time
import utils.common as common

def run_abacus(pdbfilename):
    try:
        os.mkdir("abacus_jobs")
        os.chdir("abacus_jobs")
        start_time = time.time()
        print("[INFO]: ABACUS started at %s" % (time.ctime()))
        os.system("cp ../%s ./" % (pdbfilename))
        print("[INFO]: Running ABACUS_prepare.")
        os.system("ABACUS_prepare %s" % (pdbfilename))
        print("[INFO]: Running ABACUS_S1S2.")
        os.system("ABACUS_S1S2 %s" % (pdbfilename))
        prepare_end = time.time()
        prepare_time = prepare_end - start_time
        print("[INFO]: ABACUS prepare took %f seconds." % (prepare_time))
        print("[INFO]: Running ABACUS_singleMutationScan.")

        os.system("ABACUS_singleMutationScan %s abacus_output.txt" % (pdbfilename))
        scan_end = time.time()
        scan_time = scan_end - prepare_end
        print("[INFO]: ABACUS scan took %f seconds." % (scan_time))
        os.chdir("../")
        return prepare_time, scan_time
    except FileExistsError:
        os.chdir("abacus_jobs")
        if os.path.exists("./abacus_output.txt"):
            print("[INFO]: ABACUS results found. Skipping.")
        os.chdir("../")
        return 0, 0


def runOneJob(varlist):
    def _1_2_3(x):
        d = {
            "C": "CYS",
            "D": "ASP",
            "S": "SER",
            "Q": "GLN",
            "K": "LYS",
            "I": "ILE",
            "P": "PRO",
            "T": "THR",
            "F": "PHE",
            "N": "ASN",
            "G": "GLY",
            "H": "HIS",
            "L": "LEU",
            "R": "ARG",
            "W": "TRP",
            "A": "ALA",
            "V": "VAL",
            "E": "GLU",
            "Y": "TYR",
            "M": "MET",
        }
        return d[x]

    # me varlist as foldx.runOneJob
    pdb, wild, chain, aa, resNum = varlist
    MUT = _1_2_3(aa)
    output = (
        os.popen("singleMutation %s %s %s %s" % (pdb, chain, str(resNum), MUT))
            .read()
            .split()
    )
    # print(output)
    s1 = float(output[6])
    s2 = float(output[8])
    pack = float(output[10])
    total = s1 + s2 + pack
    # self.abacus2_results["_".join([wild, str(resNum), aa])] = total
    # print(all_results)
    # A   42 GLU->TRP SAI: 0.966 S1:  1.748 S2:  0.212 PACK:  -0.009 HB:   0.000
    return "_".join([wild, str(resNum), aa]), total


def parse_abacus_out():
    try:
        os.mkdir("abacus_results")
    except FileExistsError:
        pass
    longer_names = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
    }

    with open("tempfile", "w") as tem:
        with open("abacus_jobs/abacus_output.txt") as abacusfile:
            for line in abacusfile:
                if line.startswith("site"):
                    wildAA = line.strip().split()[4]
                    wildAAnum = line.strip().split()[1]
                else:
                    tem.write(wildAA + " " + wildAAnum + " " + line)

    with open("abacus_results/All_ABACUS.score", "w+") as complete:
        complete.write(
            "#Score file formatted by GRAPE from ABACUS.\n#mutation\tscore\tstd\n"
        )
        with open("tempfile") as abacusfile:
            for line in abacusfile:
                wildAA1 = line.strip().split()[0]
                if wildAA1 in longer_names:
                    wildAAabr = longer_names[wildAA1]
                wildAAnum1 = line.strip().split()[1]
                mutAA = line.strip().split()[2]
                if mutAA in longer_names:
                    mutAAabr = longer_names[mutAA]
                sef_energy = line.strip().split()[11]
                complete.write(
                    wildAAabr
                    + "_"
                    + wildAAnum1
                    + "_"
                    + mutAAabr
                    + "\t"
                    + sef_energy
                    + "\t"
                    + str(0)
                    + "\n"
                )
            tem.close()
        complete.close()
        os.remove("tempfile")



if __name__ == "__main__":
    print("Running")
    parse_abacus_out()
