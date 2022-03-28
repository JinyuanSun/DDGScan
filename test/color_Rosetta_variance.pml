load 1pga.pdb
hide everything
show cartoon
newb=[float(i) for i in open("Rosetta_newb.txt").read().split()]
alter n. ca, b=newb.pop()
spectrum b, minimum=0, maximum=100
