# Codes for FoldX towards Stable Proteins

### See [old version](arxiv/) for previous version.

This repo will be updated within days, no more gui will be provided. However, the scoring will be more efficient and accurate since I changed `PositionScan` to `BuildModel`. In each run of `BuildModel`, 10 models will be built by default. Also in the final report, `mean ∆∆G` and `SD` will be recorded.

### Installation
First of all, please make sure you have added the foldx executable to your environment!  
The installation is simplely clone this repo and add it to PATH:
```bash
git clone https://github.com/JinyuanSun/Codes_for_FoldX.git &&
cd Codes_for_FoldX && export PATH="$(pwd):$PATH"
```
### Quickstart
You may want to try it out on a small protein like [Gb1](https://www.rcsb.org/structure/1PGA):
```bash
wget https://files.rcsb.org/download/1PGA.pdb
nohup parallelfoldscan.py -pdb 1PGA.pdb -chain A -cpu 16 -fc 1 -mode run &
```
You should expecting outputs like:  
A floder named `foldx_results` containing:
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
And another folder named `foldx_jobs` contains many subdirectories, in each subdirectory, containing raw output for every mutation built by FoldX.
### 如果你在中国大陆地区，可以使用`Gitee`:
```bash
git clone https://gitee.com/puzhunanlu30/Codes_for_FoldX.git
```
