from Bio.PDB.Polypeptide import one_to_three
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn.linear_model import LinearRegression


def run_abacus2_cmd(pdb, chain, res_id, muttype):
    cmd = f"singleMutation {pdb} {chain} {res_id} {muttype}"
    out = os.popen(cmd).read()
    lst = out.replace("\n", "").split(":")
    sai = float(lst[1].replace("S1", "").replace(" ", ""))
    s1 = float(lst[2].replace("S2", "").replace(" ", ""))
    s2 = float(lst[3].replace("PACK", "").replace(" ", ""))
    pack = float(lst[4].replace("HB", "").replace(" ", ""))
    hb = float(lst[5].replace(" ", ""))
    return [sai, s1, s2, pack, hb]


def dump_raw_data(exp_csv, pdb_dir):
    df = pd.read_csv(exp_csv)
    all_mutations = []

    for pdb, p, m in zip(df["pdb"], df["seq_index"], df["mutation"]):
        m = one_to_three(m)
        all_mutations.append(
            [os.path.join(pdb_dir, pdb + "_A.pdb"), "A", str(p + 1), m]
        )
    results = Parallel(n_jobs=20)(
        delayed(run_abacus2_cmd)(*mutation) for mutation in tqdm(all_mutations)
    )
    values = np.array(results).T
    sai, s1, s2, pack, hb = values

    df["sai"] = sai
    df["s1"] = s1
    df["s2"] = s2
    df["pack"] = pack
    df["hb"] = hb
    df.to_csv("train_abacus2.csv")
    return df


def train(abacus_df):
    coefs = []
    intercepts = []
    for i in range(10):
        regr = LinearRegression().fit(
            abacus_df[["sai", "s1", "s2", "pack", "hb"]][abacus_df["group"] != i],
            abacus_df["ddG"][abacus_df["group"] != i],
        )
        coefs.append(regr.coef_)
        intercepts.append(regr.intercept_)

    with open("~/.cache/ddgscan/abacus2_ddg_param.csv", "w+") as ofile:
        ofile.write(f"a,b,c,d,e,f\n")
        for (a, b, c, d, e), f in zip(coefs, intercepts):
            ofile.write(f"{a},{b},{c},{d},{e},{f}\n")

    def abacus2_ddg(sai, s1, s2, pack, hb):
        ddgs = []
        for (a, b, c, d, e), f in zip(coefs, intercepts):
            ddg = a * sai + b * s1 + c * s2 + d * pack + e * hb + f
            ddgs.append(ddg)
        return np.mean(ddgs), np.min(ddgs), np.std(ddgs)

    return abacus2_ddg


def get_abacus2_ddg(abacus_ddg_param_file):
    df = pd.read_csv(abacus_ddg_param_file)
    coefs = df[["a", "b", "c", "d", "e"]].values
    intercepts = df["f"]

    def abacus2_ddg(sai, s1, s2, pack, hb):
        ddgs = []
        for (a, b, c, d, e), f in zip(coefs, intercepts):
            ddg = a * sai + b * s1 + c * s2 + d * pack + e * hb + f
            ddgs.append(ddg)
        return np.mean(ddgs), np.min(ddgs), np.std(ddgs)

    return abacus2_ddg