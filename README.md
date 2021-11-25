# Towards Stable Proteins

### See [GUI](GUI/) for the GUI plugin of foldx.
GUIs only work for FoldX.

### Installation
**Warning**  
First of all, please make sure you have added the FoldX executable to your environment! Also, Rosetta is required for 
cartesian_ddg calculation or pmut_scan (a mpi build is necessary). Also, ABACUS is an outstanding software with great 
statistical energy function for protein design. Due to their licences, I cannot redistribute them here! To our glad, 
openmm is open source! So the glass is half full :).
  
To install it, clone this repo and add it to PATH:
```bash
conda install -c conda-forge openmm pdbfixer
pip install mdtraj pandas numpy joblib
git clone https://github.com/JinyuanSun/DDGScan.git &&
cd DDGScan && export PATH="$(pwd):$PATH"
```
### Usage
```
Run FoldX, Rosetta and ABACUS for in silico deep mutation scan.

arguments:
  -h, --help            show this help message and exit
  -pdb PDB, --pdb PDB   Input PDB
  -chain CHAIN, --chain CHAIN
                        Input PDB Chain to do in silico DMS
  -cpu THREADS, --threads THREADS
                        Number of threads to run FoldX, Rosetta and HHblits
  -fc FOLDX_CUTOFF, --foldx_cutoff FOLDX_CUTOFF
                        Cutoff of FoldX ddg(kcal/mol)
  -rc ROSETTA_CUTOFF, --rosetta_cutoff ROSETTA_CUTOFF
                        Cutoff of Rosetta ddg(R.E.U.)
  -ac ABACUS_CUTOFF, --ABACUS_cutoff ABACUS_CUTOFF
                        Cutoff of ABACUS SEF(A.E.U.)
  -nstruct RELAX_NUMBER, --relax_number RELAX_NUMBER
                        Number of how many relaxed structure
  -nruns NUMOFRUNS, --numofruns NUMOFRUNS
                        Number of runs in FoldX BuildModel
  -sl SOFTLIST, --softlist SOFTLIST
                        List of Software to be used in ddg calculations
  -mode MODE, --mode MODE
                        Run, Rerun or analysis
  -preset PRESET, --preset PRESET
                        Fast or Slow
  -md MOLECULAR_DYNAMICS, --molecular_dynamics MOLECULAR_DYNAMICS
                        Run 1ns molecular dynamics simulations for each mutation using openmm.
  -platform PLATFORM, --platform PLATFORM
                        CUDA or CPU

```


### Quickstart
You may want to try it out on a small protein like [Gb1](https://www.rcsb.org/structure/1PGA):  
I will recommend to use the `fast` mode with `-md` turned on `True` and using `CUDA` to accelerate molecular dynamics simulations.
```bash
wget https://files.rcsb.org/download/1PGA.pdb
parallel_single_scan.py -pdb 1PGA.pdb -chain A -sl "FoldX,Rosetta,ABACUS" -mode run -cpu 40 -preset fast -md True -platform CUDA
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
every mutation built by FoldX. Of course there will be folder start with rosetta or abacus, depending on your choice!  
If `-md` was turned on, all produced snapshots can by found in `selectpdb` with `afterMD` as suffix in the name of pdbfiles.
### Inspect structures
Using `scripts/inspectmutation.py` to inspect mutations in pymol:
```bash
pymol inspectmutation.py $Wildtype_structure $Mutation_structure $Mutation_position $Chain
```
### Citations
If you find the models useful in your research, we ask that you cite the relevant paper:

```bibtex

@article{cui20201cascatal,
    author = {Cui, Yinglu and Chen, Yanchun and Liu, Xinyue and Dong, Saijun and Tian, Yu’e and Qiao, Yuxin and Mitra, Ruchira and Han, Jing and Li, Chunli and Han, Xu and Liu, Weidong and Chen, Quan and Wei, Wangqing and Wang, Xin and Du, Wenbin and Tang, Shuangyan and Xiang, Hua and Liu, Haiyan and Liang, Yong and Houk, Kendall N. and Wu, Bian},
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
	Author = {Sun, Jinyuan and Cui, Yinglu and Wu, Bian},
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


### 如果你在中国大陆地区，可以使用`Gitee`:
```bash
git clone https://gitee.com/puzhunanlu30/Codes_for_FoldX.git
```
