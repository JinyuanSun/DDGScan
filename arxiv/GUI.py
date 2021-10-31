import os
import ssl
import urllib
import urllib.request
from tkinter import *
from tkinter import messagebox

ssl._create_default_https_context = ssl._create_unverified_context
'''
LARGE_FONT= ("Verdana", 12)
NORM_FONT = ("Helvetica", 10)
SMALL_FONT = ("Helvetica", 8)
'''
os.system("mkdir EasyFoldx")
sourcePath = "EasyFoldx"
os.chdir(sourcePath)

def mutation():
    inputname = pdb_text.get()
    chainid = pdb_text.get()[5]
    pdbname = inputname.split(".")[0] + "_Repair.pdb"
    mut = mutation_text.get()
    l = mut[0] + chainid + mut[1:]
    print(l)
    o = "command=PositionScan\npdb=" + pdbname + "\npositions=" + l
    of_name = "MT_" + mut + ".cfg"
    of = open(of_name, "w+")
    of = open(of_name, "a+")
    print(o, file=of)
    #cmd = "./foldx --config " + of_name
    #os.system("./foldx -f " + of_name)
    os.system("nohup ./foldx -f " + of_name+" &")


def analyze():
    inputname = pdb_text.get()
    filename = "PS_"+inputname+"_Repair_scanning_output.txt"
    '''
    file = open(filename)
    l = ''
    for x in file:
        l = l +x+"\n"
    mut_out.insert(END,l)
    '''
    data = []
    with open(filename) as f:
        for line in f:
            data += line.split()

#Create your listbox here.
    for i in range(len(data)):
        mut_out.insert(i+1, data[i])

    #print("./foldx -f " + of_name)
#

def PS():
    inputname = pdb_text.get()
    pdbname = inputname.split(".")[0] + "_Repair.pdb"
    print(pdbname)
    # pdbname = args.pdbfile
    nt = int(nt_text.get())

    # nt = 20
    # pdbname="6QG9_A_Repair.pdb"
    import os
    try:
        file = open("SO_" + pdbname.replace("pdb", "fxout"), "r")
    except FileNotFoundError:
        os.system("./foldx --command=SequenceOnly --pdb=" + pdbname)
        file = open("SO_" + pdbname.replace("pdb", "fxout"), "r")
    lst = []
    for line in file:
        l = line.replace("\n", "").split("\t")
        if len(l) > 3:
            lst.append(l[3] + "a")

    t = len(lst) // (nt - 1)
    n = 0
    for i in range(0, len(lst), t):
        b = lst[i:i + t]
        l = ""
        for x in b:
            l = l + x + ","
        n = n + 1
        o = "command=PositionScan\npdb=" + pdbname + "\npositions=" + l
        of_name = "PS_" + str(n) + ".cfg"
        of = open(of_name, "w+")
        of = open(of_name, "a+")
        print(o[:-1], file=of)
        os.system("nohup ./foldx -f " + of_name + " &")


def populate_list():
    print("populate")


def getpdb():
    if len(pdb_text.get()) < 5:
        pdbid = pdb_text.get()
        pdb = pdbid + ".pdb"
        urllib.request.urlretrieve('https://files.rcsb.org/download/' + pdb, pdb)
        cwd = os.getcwd() + "/" + pdb
        messagebox.showinfo("Done", "The pdb file " + pdbid + " has been downloaded to " + cwd)
        oname = pdbid
    else:
        pdbid = pdb_text.get()[0:4]
        pdb = pdbid + ".pdb"
        urllib.request.urlretrieve('https://files.rcsb.org/download/' + pdb, pdb)
        pdbfile = open(pdb)
        oname = pdb_text.get() + ".pdb"
        opdb = open(oname, "w+")
        opdb = open(oname, "a+")
        chainid = pdb_text.get()[5]
        for line in pdbfile:
            if line.startswith("ATOM"):
                l = line.split()
                if l[4] == chainid:
                    print(line, end="", file=opdb)
        cwd = os.getcwd() + "/" + oname
        messagebox.showinfo("Done", "The pdb file " + oname + " has been downloaded to " + cwd)
    return oname
    # os.system("")


def repair_pdb():
    inputname = pdb_text.get()+".pdb"
    pdbname = inputname.split(".")[0] + "_Repair.pdb"
    os.system("./foldx --command=RepairPDB --pdb=" + inputname)
    cwd = os.getcwd() + "/" + pdbname
    messagebox.showinfo("Done", "The Repaired pdb file " + pdbname + " has been saved to " + cwd)


# create window object
app = Tk()

# pdb
pdb_text = StringVar()
pdb_lable = Label(app, text='PDB ID', font=('bold', 14), pady=20)
pdb_lable.grid(row=0, column=0, sticky=W)
pdb_entry = Entry(app, textvariable=pdb_text)
pdb_entry.grid(row=0, column=1)

# mutation
mutation_text = StringVar()
mutation_lable = Label(app, text='Mutation', font=('bold', 14))
mutation_lable.grid(row=2, column=0, sticky=W)
mutation_entry = Entry(app, textvariable=mutation_text)
mutation_entry.grid(row=2, column=1)

# chain
inpdb_text = StringVar()
inpdb_lable = Label(app, text='Input PDB', font=('bold', 14))
inpdb_lable.grid(row=1, column=0, sticky=W)
inpdb_entry = Entry(app, textvariable=pdb_text)
inpdb_entry.grid(row=1, column=1)

# PS
nt_text = StringVar()
nt_lable = Label(app, text='Number of Threads', font=('bold', 14))
nt_lable.grid(row=3, column=0, sticky=W)
nt_entry = Entry(app, textvariable=nt_text)
nt_entry.grid(row=3, column=1)


# mutation out
mut_out = Listbox(app, height=8, width=50, border=0)
mut_out.grid(row=4, column=0, columnspan=3, rowspan=6, pady=20, padx=20)

# creat scrollbar
scrollbar = Scrollbar(app)
scrollbar.grid(row=4, column=3)


# Buttons
fetch_btn = Button(app, text='Fetch', width=12, command=getpdb)
fetch_btn.grid(row=0, column=2, )

repair_btn = Button(app, text='Foldx Repair', width=12, command=repair_pdb)
repair_btn.grid(row=1, column=2)

Position_Scan_btn = Button(app, text='Foldx Position Scan', width=15, command=PS)
Position_Scan_btn.grid(row=3, column=2)

mutation_btn = Button(app, text='Foldx Mutation', width=12, command=mutation)
mutation_btn.grid(row=2, column=2)

analyze_btn = Button(app, text='Show', width=12, command=analyze)
analyze_btn.grid(row=5, column=0)


def pwd():
    print(os.getcwd())

pwd_btn = Button(app, text='pwd', width=12, command=pwd)
pwd_btn.grid(row=2, column=3)
# set scroll to listbox
# mut_out.configure(yscrollcommand=scrollbar.set)
# scrollbar.configure(command=mut_out.yview)


app.title('Easy FoldX')
app.geometry('700x500')

# populate
populate_list()

# start program
app.mainloop()
