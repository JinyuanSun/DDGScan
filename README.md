# Codes_for_FoldX
some useful codes for [foldx](http://foldxsuite.crg.eu/products#foldx) 
# [pdb2seq](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/pdb2seq.py)
convert a pdb file to a fasta  
usage: python pdb2seq.py name.pdb
# [prepare4scan](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/prepare4scan.py)
to build the config file for foldx to preform a saturated mutation scan  
usage: python prepare4scan.py -i seq.fasta -p name_Repair.pdb
# prepare4buildmodel
prepare the individual list for BuildModel in FoldX   
usage:python prepare4buildmodel.py -i PS_pdbid_scanning_output.txt -l a_number
