import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='prepare for BuildModel in FoldX')
parser.add_argument("-i", '--input', help="input a txt file, a PS_pdbid_file.txt mostly")
parser.add_argument("-l", '--lines', help='number of lines you want')

args = parser.parse_args()

inf = args.input
lines = args.lines


def _3_2_1(x):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
         'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    y = d[x]
    return y


def multi_sub(self,string,p,c):
    new = []
    for s in string:
        new.append(s)
    for index,point in enumerate(p):
        new[point] = c[index]
    return ''.join(new)


data = pd.read_table(inf, error_bad_lines=False, header=None)
sort_data = data.sort_values(by=[1])
part_data = sort_data.head(int(lines))

ofile = open('individual_list.txt', "w")
_list = np.array(part_data[0]).tolist()
for i in _list:
    _3 = i[0:3]
    mut = _3_2_1(_3) + i[3:len(i)] + ';'
    print(mut,file=ofile)
