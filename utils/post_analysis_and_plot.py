#!/usr/bin/env python     
# -*- coding: utf-8 -*-
# @Author  : Jinyuan Sun
# @Time    : 2022/3/27 4:20 PM
# @File    : post_analysis_and_plot.py
# @annotation    :

# TODO: venn, heatmap, position_avg, write scaled avg ddg to b factor

import argparse
import distutils.dir_util

import logomaker as lm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from venn import venn


def score_file2array(all_score_file):
    """
    param: all_score_file (FoldX, Rosetta and ABACUS)
    return: heatmap array
    """
    global method
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    aa2num = {}
    energy_dict = {}
    resnum_list = []
    resnum2index = {}
    wild_type_dict = {}

    score_file = open(all_score_file, 'r')

    for i, aa in enumerate(alphabet):
        aa2num[aa] = i

    for line in score_file:
        if line[0] != '#':
            mutation, ddg, std = line.replace("\n", "").split("\t")
            energy_dict[mutation] = (float(ddg), float(std))
            resnum = int(mutation.split("_")[1])
            resnum_list.append(resnum)
            wild_type_dict[resnum] = mutation.split("_")[0]
        if line.startswith("#Score"):
            method = line.replace("\n", "").split(' ')[-1][:-1]
    score_file.close()

    total_length = len(set(resnum_list))
    for i, resnum in enumerate(np.sort(list(set(resnum_list)))):
        resnum2index[resnum] = i

    ddg_array = np.zeros((20, total_length))
    std_array = np.zeros((20, total_length))
    for key, value in energy_dict.items():
        i = aa2num[key[-1]]
        j = resnum2index[int(key.split("_")[1])]
        ddg_array[i][j] = value[0]
        std_array[i][j] = value[1]
    return ddg_array, std_array, list(resnum2index.keys()), list(resnum2index.values()), method, wild_type_dict


def energy2logo(arr, resnum_list, index_list, method):
    """
    param: energy_array: (20,)
    p_mut/p_ref = exp( beta * -ddG)
    """
    plt.figure(figsize=(arr.shape[1] / 5, arr.shape[0] / 5), dpi=150)
    p_ratio = np.exp(-arr.T / 0.6)
    z = np.sum(p_ratio, axis=0)
    p_ref = 1 / z
    #
    pesudo_pssm = p_ratio * p_ref
    data = pd.DataFrame(index=index_list,
                        data=pesudo_pssm / np.sum(pesudo_pssm, axis=1).reshape(-1, 1),
                        columns=list("ARNDCQEGHILKMFPSTWYV"))
    res_logo = lm.Logo(data,
                       fade_probabilities=True,
                       stack_order='small_on_top',
                       font_name='Arial Rounded MT Bold')
    res_logo.ax.set_xticks(np.array(index_list))
    res_logo.ax.set_title(f"{method}")
    res_logo.ax.set_xticklabels(resnum_list, rotation='vertical')
    plt.savefig(f'plots/{method}_Logo.png', bbox_inches='tight')
    # p_else.shape
    # np.sum(p_else, axis=0) + p_ref
    # return data


def heatmap(arr, resnum_list, index_list, method):
    """
    param: arr: array shape (20, L)
    """

    plt.figure(figsize=(arr.shape[1] / 5, arr.shape[0] / 5), dpi=150)
    ax = sns.heatmap(arr, cmap='Blues',
                     cbar_kws={"orientation": "vertical", "pad": 0.01})
    ax.set_xticks(np.array(index_list) + .5)
    ax.set_xticklabels(resnum_list, rotation='vertical')
    ax.set_yticks(np.array(range(0, 20)) + .5)
    ax.set_yticklabels(list("ARNDCQEGHILKMFPSTWYV"))
    ax.set_title(f"{method}")
    plt.savefig(f'plots/{method}_ddg.png', bbox_inches='tight')


def position_avg(arr, resnum_list, index_list, method):
    """
    param: arr: array shape (20, L)
    """
    plt.figure(figsize=(arr.shape[1] / 5, arr.shape[0] / 5), dpi=150)
    arr = arr.reshape(-1, )
    position = np.repeat(resnum_list, 20)
    indexes = np.repeat(index_list, 20)
    data = pd.DataFrame(index=indexes, data={'position': position, 'ddg': arr})
    data = data.dropna()
    ax = sns.boxplot(x=data["position"], y=data["ddg"], palette="Greys")
    ax.set_xticks(np.array(index_list))
    ax.set_xticklabels(resnum_list, rotation='vertical')
    ax.set_title(f"{method}")
    plt.savefig(f'plots/{method}_boxplot.png', bbox_inches='tight')


def kde_plot(arr, method):
    """
    param: arr: array shape (20, L)
    """
    plt.figure(figsize=(4, 4), dpi=150)
    arr = arr.reshape(-1, )

    ax = sns.kdeplot(arr)
    ax.set_title(f"{method}")
    plt.savefig(f'plots/{method}_kde.png', bbox_inches='tight')


def posistion_variance(arr, resnum_list, index_list, method):
    """
    param: arr: array shape (20, L)
    """
    plt.figure(figsize=(arr.shape[1] / 5, arr.shape[0] / 5), dpi=150)
    variance = np.std(arr, axis=0)

    ax = sns.lineplot(x=index_list, y=variance, label="Variances", marker='o')
    ax.set_xticks(np.array(index_list))
    # ax.legend()
    ax.set_xticklabels(resnum_list, rotation='vertical')
    ax.set_xlim(index_list[0], index_list[-1])
    ax.set_title(f"{method}")
    plt.savefig(f'plots/{method}_variance.png', bbox_inches='tight')


def venn_plot(engine_cutoff_dict):
    """
    param: engine_cutoff_dict: {"rosetta":1.5,...}
    """
    alphabet = "ARNDCQEGHILKMFPSTWYV"
    aa2num = {}
    engine_selects = {}

    for i, aa in enumerate(alphabet):
        aa2num[aa] = i

    for engine, cutoff in engine_cutoff_dict.items():
        selects = []
        score_file = open(f"{engine.lower()}_results/All_{engine}.score", 'r')
        for line in score_file:
            if line[0] != '#':
                mutation, ddg, std = line.replace("\n", "").split("\t")
                if float(ddg) <= cutoff:  # cutoff should be converted to negative
                    selects.append(mutation)
        score_file.close()
        engine_selects[engine] = set(selects)

    fig, ax = plt.subplots(figsize=(5, 6), dpi=150)

    venn(engine_selects, ax=ax, cmap='Set3')
    ax.legend(list(engine_selects.keys()), ncol=3, loc=(.1, 0.9))
    #            #            'MSAddg %d/%d'%(len(set(msa_selected)), msa_predicted),
    #            'ABACUS %d/%d' % (len(set(abacus_selected)), abacus_predicted),
    #            'FoldX %d/%d' % (len(set(foldx_selected)), foldx_predicted)], loc=(.15, 0.9))
    plt.savefig('plots/venn.png', bbox_inches='tight')


def residue_bar(engine_cutoff_dict, residue_position):
    """
    param: engine_cutoff_dict: {"rosetta":-1.5,...}
    """

    global wild_type
    engine_selects = {'Mutation': [],
                      'Predicted score': [],
                      'Engine': []}

    for engine, cutoff in engine_cutoff_dict.items():
        score_file = open(f"{engine.lower()}_results/All_{engine}.score", 'r')
        for line in score_file:
            if line.startswith("#Score"):
                method = line.replace("\n", "").split(' ')[-1][:-1]
            if line[0] != '#':
                mutation, ddg, std = line.replace("\n", "").split("\t")
                if str(residue_position) == mutation.split("_")[1]:
                    wild_type, position, mut_type = mutation.split("_")
                    engine_selects['Mutation'].append(mut_type)
                    engine_selects['Predicted score'].append(float(ddg))
                    engine_selects['Engine'].append(method)
        score_file.close()

    data = pd.DataFrame(data=engine_selects)
    fig, ax = plt.subplots(figsize=(8, 5), dpi=150)
    ax = sns.barplot(x='Mutation', y='Predicted score', hue="Engine", data=data, palette="Blues")
    ax.set_title(f'{wild_type}{residue_position}_histogram')
    ax.set_ylabel(f'Scores')
    ax.set_xlabel('Mutations')
    plt.savefig(f'plots/{wild_type}{residue_position}_histogram.png', bbox_inches='tight')


def write_variance2ca(arr, method, pdb):
    variance = np.mean(arr, axis=0)
    variance = (100 * (variance - np.min(variance)) / np.ptp(variance)).astype(int)
    with open(f"{method}_newb.txt", 'w+') as outfile:
        for x in variance:
            outfile.write(f'{x} ')
        outfile.close()
    with open(f"color_{method}_variance.pml", "w+") as cmdfile:
        cmdfile.write(f"load {pdb}\n")
        cmdfile.write("hide everything\nshow cartoon\n")
        cmdfile.write(f'newb=[float(i) for i in open("{method}_newb.txt").read().split()]\n')
        cmdfile.write("alter n. ca, b=newb.pop()\n")
        cmdfile.write(f"spectrum b, minimum={min(variance)}, maximum={max(variance)}\n")
        # cmdfile.write("as cartoon\ncartoon putty\nset cartoon_putty_radius, 1")


def main(args):
    residue_position = args.residue_position  # = 26
    distutils.dir_util.mkpath('plots')
    engine_cutoff_dict = args.engine_cutoff_dict  # = {'rosetta': 0, 'FoldX': 0, 'ABACUS': 0}
    pdb = args.pdb  # '1pga.pdb'
    plot_type = args.plot_type
    for engine in engine_cutoff_dict:
        all_score_file = f"{engine.lower()}_results/All_{engine}.score"
        arr, std, resnum_list, index_list, method, wild_type_dict = score_file2array(all_score_file)
        if plot_type == 'all':
            write_variance2ca(arr, method, pdb)
            heatmap(arr, resnum_list, index_list, method)
            position_avg(arr, resnum_list, index_list, method)
            posistion_variance(arr, resnum_list, index_list, method)
            kde_plot(arr, method)
            energy2logo(arr, resnum_list, index_list, method)
        if 'venn' in plot_type:
            venn_plot(engine_cutoff_dict)
        if 'residue_bar' in plot_type and residue_position:
            residue_bar(engine_cutoff_dict, residue_position)
        if 'heatmap' in plot_type:
            heatmap(arr, resnum_list, index_list, method)
        if 'position_avg_boxplot' in plot_type:
            position_avg(arr, resnum_list, index_list, method)
        if 'variance_lineplot' in plot_type:
            posistion_variance(arr, resnum_list, index_list, method)
        if 'kde_plot' in plot_type:
            kde_plot(arr, method)
        if 'residue_logo' in plot_type:
            energy2logo(arr, resnum_list, index_list, method)


def get_args():
    parser = argparse.ArgumentParser(
        description="analysis and plot results from grape_phaseI result"
    )
    parser.add_argument("pdb", help="your target pdb file")
    parser.add_argument("results_dir", help="directory of results of grape_phase_I or list_distribute")
    parser.add_argument("--residue_position", help="residue position, if you asked for a barplot at residue level")
    parser.add_argument("--plot_type",
                        help="plots you want to make",
                        choices=["all",
                                 "venn",
                                 'residue_bar',
                                 'heatmap',
                                 'position_avg_boxplot',
                                 'variance_lineplot',
                                 'kde_plot',
                                 'residue_logo'
                                 ]
                        )

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    main(args)

