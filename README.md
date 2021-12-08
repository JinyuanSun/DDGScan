# Towards Stable Proteins

- [Towards Stable Proteins](#towards-stable-proteins)
    + [The GUI plugin for FoldX.](#the-gui-plugin-for-foldx)
    + [Installation](#installation)
    + [Usage](#usage)
    + [QuickStart](#quickstart)
    + [Inspect structures](#inspect-structures)
    + [Citations](#citations)
    + [Others](#others)

**I am testing this repo with some different input structures, if you encountered any failure please post a issue.** 

### The GUI plugin for FoldX.
[GUI](GUI/) only work for FoldX.

### Installation
 
First of all, please make sure you have added the **FoldX** executable to your environment! Secondly, **Rosetta** 
(a mpi build is necessary) is 
required for cartesian_ddg (`-mode slow`) calculation or pmut_scan(`-mode fast`). 
Also, **ABACUS** is an outstanding software with great statistical energy function for protein design. 
Structures downloaded from RCSB could be erroneous. One of the biggest problems that will directly affect energy calculation is breaks in chains. 
Here I implemented a loop closure module using **modeller**, a great software with a very long history, as backend.   
Due to their licenses, I cannot redistribute them here :worried: !  
To our glad, **openmm** is open source! So the glass is half full :smiley: .
  
To install it, clone this repo and add it to PATH:
```bash
conda install -c conda-forge openmm pdbfixer
pip install mdtraj pandas numpy joblib
git clone https://github.com/JinyuanSun/DDGScan.git &&
cd DDGScan && export PATH="$(pwd):$PATH"
```
### Usage
I provide many options for users especially those know what they want. I really tried to make this package light and also 
be well functional. Here are some quick walk-through. `pdb` and `chain` are positional but really you need to set 
`-E` according to the software you have in your OS. `-seq` are strongly recommended to be set by the user. 
Also, I highly recommend adding the `-MD` flag and using `-P CUDA` if a good gpu is available (better
 than RTX2060 well be much faster than 48 core cpu). Also, I did not test how much precision dropped to use the `-S fast` 
 preset, but I do know it can be faster in about two orders of magnitude.  
 If using `-fill` flag, input structure will be automatically fixed using information from SEQRES record in native PDB 
 downloaded from RCSB using modeller. Model with lowest `molpdf` energy will be subjected to following step.
```
Run FoldX, Rosetta and ABACUS for in silico deep mutation scan.

positional arguments:
  pdb                   Input PDB
  chain                 Input PDB Chain to do in silico DMS

optional arguments:
  -h, --help            show this help message and exit
  -fill, --fill_break_in_pdb
                        Use modeller to fill missing residues in your pdb file. Use this option with caution!
  -seq SEQUENCE, --sequence SEQUENCE
                        The exact sequence of protein you want to design. All mutation will be named according to this sequence.
  -T THREADS, --threads THREADS
                        Number of threads to run FoldX, Rosetta
  -fc FOLDX_CUTOFF, --foldx_cutoff FOLDX_CUTOFF
                        Cutoff of FoldX ddg(kcal/mol)
  -rc ROSETTA_CUTOFF, --rosetta_cutoff ROSETTA_CUTOFF
                        Cutoff of Rosetta ddg(R.E.U.)
  -ac ABACUS_CUTOFF, --abacus_cutoff ABACUS_CUTOFF
                        Cutoff of ABACUS SEF(A.E.U.)
  -nstruct RELAX_NUMBER, --relax_number RELAX_NUMBER
                        Number of how many relaxed structure
  -nruns NUMOFRUNS, --numofruns NUMOFRUNS
                        Number of runs in FoldX BuildModel
  -E {abacus,foldx,rosetta} [{abacus,foldx,rosetta} ...], --engine {abacus,foldx,rosetta} [{abacus,foldx,rosetta} ...]
  -M {run,rerun,analysis,test}, --mode {run,rerun,analysis,test}
                        Run, Rerun or analysis
  -S {fast,slow}, --preset {fast,slow}
                        Fast or Slow
  -MD, --molecular_dynamics
                        Run 1ns molecular dynamics simulations for each mutation using openmm.
  -P {CUDA,CPU}, --platform {CUDA,CPU}
                        CUDA or CPU
```


### QuickStart
You may want to try it out on a small protein like [Gb1](https://www.rcsb.org/structure/1PGA):  
I will recommend using the `-S fast` with `-MD` flag, and using `CUDA` to accelerate molecular dynamics simulations. 
This is a very good crystal structure solved by X-ray, so I did not pass any value about fixing the PDB file!  
Using `-S slow` to get more accuracy!
```bash
wget https://files.rcsb.org/download/1PGA.pdb
grape-fast.py 1PGA.pdb A -E "FoldX,Rosetta,ABACUS" -M run -T 40 -S slow -MD -P CUDA
```
You should expecting outputs like:  
A folder named `foldx_results` containing:
```
All_FoldX.score
MutationsEnergies_BestPerPositionBelowCutOff_SortedByEnergy.tab
MutationsEnergies_BelowCutOff.tab
MutationsEnergies_BestPerPosition_SortedByEnergy.tab
MutationsEnergies_BelowCutOff_SortedByEnergy.tab
MutationsEnergies_CompleteList.tab
MutationsEnergies_BestPerPosition.tab
MutationsEnergies_CompleteList_SortedByEnergy.tab
MutationsEnergies_BestPerPositionBelowCutOff.tab
```
And another folder named `foldx_jobs` contains many subdirectories, in each subdirectory, containing raw output for 
every mutation built by FoldX. Of course, there will be directories start with rosetta or abacus, depending on your choice!  
If `-md` was turned on, all produced snapshots can be found in `selectpdb` with `afterMD` as a suffix in the name of PDB files.
### Inspect structures
Using `scripts/inspectmutation.py` to inspect mutations in pymol:
```bash
pymol inspectmutation.py $Wildtype_structure $Mutation_structure $Mutation_position $Chain
```
### Citations
If you find the models useful in your research, we ask that you cite the relevant paper:

```bibtex
@article {cui2021cascatal,
    author = {Cui YL, Chen YC, Liu XY, Dong SJ, Tian YE, Qiao YX, Mitra R, Han J, Li CL, Han X, Liu WD, Chen Q, Wei WQ, Wang X, Du, Tang SY, Xiang H, Liu HY, Liang Y, Houk KN, Wu B},
    title = {Computational Redesign of a PETase for Plastic Biodegradation under Ambient Condition by the GRAPE Strategy},
    journal = {ACS Catalysis},
    volume = {11},
    number = {3},
    pages = {1340-1350},
    year = {2021},
    doi = {10.1021/acscatal.0c05126}
}
@article {sun2021mie,
    Title = {GRAPE, a greedy accumulated strategy for computational protein engineering},
    Author = {Sun JY and Cui YL and Wu B},
    DOI = {10.1016/bs.mie.2020.12.026},
    Volume = {648},
    Year = {2021},
    Journal = {Methods in enzymology},
    ISSN = {0076-6879},
    Pages = {207—230},
    URL = {https://doi.org/10.1016/bs.mie.2020.12.026}
}
@article {cui2022,
    Title = {GRAPE-web: A web server for automated design of thermostable proteins. (in prep.)},
    Author = {Cui, Yinglu and Wu, Bian}
}
```


### Others
如果你在中国大陆地区，可以使用
```bash
git clone https://gitee.com/puzhunanlu30/Codes_for_FoldX.git
```
or try this:
```bash
git clone https://github.com.cnpmjs.org/JinyuanSun/DDGScan.git
```

### Develop Information
2021.10: Restart this project.  
2021.11: Added `openmm` for MDs.  
2021.12: Added `modeller` for loop modelling and args was rewritten.  
Developed this in every day 20:00 - 02:00 :cat: . Continuing...
