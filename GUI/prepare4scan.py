import argparse

parser = argparse.ArgumentParser(description='prepare for PositionScan in FoldX')
parser.add_argument("-i", help="input a fasta file")
parser.add_argument("-p", help="the a pdb file corresponding to the fasta file, a repaired pdbfile mostly")

args = parser.parse_args()

inf = open(args.i)
pdbid = args.p


def readseq(fasta):
    seq = ''
    for line in fasta:
        line = line.strip('\n')
        if not line.startswith('>'):
            seq += line
    return seq


ofile = open('config_scan.cfg', "w")
print('command=PositionScan' + '\n' + 'pdb=' + pdbid + '\n' + 'positions=', end='', file=ofile)

seq = readseq(inf)
for i in range(len(seq)):
    if i < len(seq) - 1:
        new_word = seq[i] + 'A' + str(i + 1) + 'a,'
        print(new_word, end='', file=ofile)
    else:
        new_word = seq[i] + 'A' + str(i + 1) + 'a'
        print(new_word, end='', file=ofile)
