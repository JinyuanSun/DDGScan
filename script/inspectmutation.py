import pymol
import os
import sys

wt = sys.argv[2]
mut = sys.argv[3]
pos = sys.argv[4]
chain = sys.argv[5]
wt_name = wt.replace(".pdb","")
mut_name = mut.replace(".pdb","")
cmd.load(wt)
cmd.load(mut)
cmd.remove("resn HOH")
cmd.show("sticks")
cmd.distance("hbonds",mut_name,mut_name,"3.2","2")
cmd.distance("hbonds",wt_name,wt_name,"3.2","2")
cmd.hide("labels","all")
cmd.util.cbag(wt_name)
cmd.util.cbab(mut_name)
cmd.util.cbac("/"+wt_name+"//"+chain+"/"+pos+"/")
cmd.util.cbap("/"+mut_name+"//"+chain+"/"+pos+"/")
cmd.select("reg",wt_name+"//"+chain+"/"+pos+"/ expand 8")
cmd.zoom("reg")
cmd.delete("reg")
cmd.space("CMYK")
cmd.bg_color("white")
