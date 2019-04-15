import sys


def _get_sequence(pdb_file):
    '''Extract the sequence from a PDB file'''
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import three_to_one

    s = PDBParser().get_structure("id", pdb_file)

    sequence = "".join([three_to_one(x.get_resname()) for x in s.get_residues()])

    return sequence


pdbfile = sys.argv[1]
seqID = pdbfile[0:len(pdbfile) - 4]
fasfile = open(seqID + '.fasta', "w")
print('>' + seqID, file=fasfile)
print(_get_sequence(pdbfile), file=fasfile)
