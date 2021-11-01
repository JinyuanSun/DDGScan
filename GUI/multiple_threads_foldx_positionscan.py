import argparse
parser = argparse.ArgumentParser(description=
                                 'To run Foldx PositionScan with multiple threads, make sure'+
                                 ' that you have the foldx and your pdb in the same floder')
parser.add_argument("-s", '--pdbfile', help="The pdb file, the repiared one")
parser.add_argument("-nt", '--number_threads', help="How many threads to run the Foldx")
args = parser.parse_args()

pdbname = args.pdbfile
nt = args.number_threads

#nt = 20
#pdbname="6QG9_A_Repair.pdb"
import os
try:
    file = open("SO_"+pdbname.replace("pdb","fxout"),"r")
except FileNotFoundError:
    os.system("./foldx --command=SequenceOnly --pdb="+pdbname)
    file = open("SO_"+pdbname.replace("pdb","fxout"),"r")
lst = []
for line in file:
    l = line.replace("\n","").split("\t")
    if len(l)>3:
        lst.append(l[3]+"a")
        
t = len(lst)//(int(nt)-1)
n = 0
for i in range(0,len(lst),t):
    b=lst[i:i+t]
    l = ""
    for x in b:
        l = l+x+","
    n = n + 1
    o = "command=PositionScan\npdb="+pdbname+"\npositions="+l
    of_name = "PS_"+str(n)+".cfg"
    of = open(of_name,"w+")
    of = open(of_name,"a+")
    print(o[:-1],file=of)
    os.system("nohup ./foldx -f "+of_name+" &")
