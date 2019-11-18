# Codes_for_FoldX
some useful codes for [foldx](http://foldxsuite.crg.eu/products#foldx) 
# [pdb2seq](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/pdb2seq.py)
convert a pdb file to a fasta  
usage: python pdb2seq.py name.pdb
# [prepare4scan](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/prepare4scan.py)
to build the config file for foldx to preform a saturated mutation scan  
usage: python prepare4scan.py -i seq.fasta -p name_Repair.pdb
# [prepare4buildmodel](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/prepare4buildmodel.py)
prepare the individual list for [BuildModel](http://foldxsuite.crg.eu/command/BuildModel) in FoldX   
usage:python prepare4buildmodel.py -i PS_pdbid_scanning_output.txt -l a_number
# [multiple_threads_positionscan](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/multiple_threads_foldx_positionscan.py)
To run the positionscan in a faster way with more threads   
usage:python multiple_threads_foldx_positioncsan.py -s pdb -nt number_of_threads

# [GUI script](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/GUI.py)
For Biologists who has difficulty in command-line tools.  
This script can be simply builded by: `pyinstaller GUI.py`  

# Guide for the most people  
1. Download the GUI.py  
2. make sure you have tkinter installed, which can be install by the following command:  
  - LINUX: `sudo apt-get install python3-tk`
  - MAC: Most already installed and you can refer to this link: https://discourse.brew.sh/t/python3-installation-with-tkinter/3636  
3. run `python3 GUI.py`  
4. Usage  
  - After the PDB ID, you can enter the four-digit PDB code, or you can enter six digits: 3wzl_A represents the A chain of  3wzl, and the acquisition is directly downloaded from the RCSB database.  
  - At present, the PDB ID and the input PDB are synchronized and can be modified by user. Foldx Repair is a preparation step for calculating the mutation ddG, and the calculation of A chain of 3wzl takes about several minutes.  
  - Mutation currently only supports single point mutations (eg. M1A, mutate the first Met to Ala)  
  - The <"Foldx Position Scan"> is currently not recommended. This is a saturated mutation scanner, the input is the number of threads, because foldx does not support MPI, so this is a manual multi-thread. See [multiple_threads_positionscan](https://github.com/github.com/JinyuanSun/Codes_for_FoldX/blob/master/multiple_threads_foldx_positionscan.py) for details. 
If you have any questions, please contact 1650949260@qq.com  
# GUI版本foldx使用教程（目前支持MAC,LINUX)  
# 因为foldx本身是需要许可证的，所以我不能在这里提供下载，以下的说明是针对我所在的组里面的成员的。
1.安装  
  - 打开终端  
  - 输入： mkdir EasyFoldx  
  - 将压缩包中的foldx和rotabase.txt拷贝入EasyFoldx  
  - 在dist文件夹内双击part_manger  
2. 使用  
  - PDB ID 后可输入四位PDB code，也可以输入六位：3wzl_A即代表3wzl的A链，fetch是从RCSB数据库直接下载。  
  - 目前PDB ID和Input PDB是同步的，可自行修改，Foldx Repair是计算突变ddG的准备步骤，3wzl的A链大概需要几分钟。  
  - Mutation目前仅支持单点突变，（如M1A，即为第一位的Met突变为Ala）  
  - Number of Threads目前不建议使用，这是一个饱和突变扫描程序，输入为线程数，因为foldx不支持MPI，所以这是一个手动多线程。具体内容可参见[multiple_threads_positionscan](https://github.com/JinyuanSun/Codes_for_FoldX/blob/master/multiple_threads_foldx_positionscan.py)。  
如有任何问题，请联系1650949260@qq.com
# GUI version foldx tutorial (currently supports MAC and LINUX only)  
# Beacues Foldx requires lisence, I can't provide it here, the following guide is for members in the lab I'm working with.  
1. Installation  
  - Open terminal  
  - Enter: mkdir EasyFoldx  
  - Copy foldx and rotabase.txt from the archive to EasyFoldx  
  - Double-click part_manger in the dist folder  
2. Usage  
  - After the PDB ID, you can enter the four-digit PDB code, or you can enter six digits: 3wzl_A represents the A chain of  3wzl, and the acquisition is directly downloaded from the RCSB database.  
  - At present, the PDB ID and the input PDB are synchronized and can be modified by user. Foldx Repair is a preparation step for calculating the mutation ddG, and the calculation of A chain of 3wzl takes about several minutes.  
  - Mutation currently only supports single point mutations (eg. M1A, mutate the first Met to Ala)  
  - The <"Foldx Position Scan"> is currently not recommended. This is a saturated mutation scanner, the input is the number of threads, because foldx does not support MPI, so this is a manual multi-thread. See [multiple_threads_positionscan](https://github.com/github.com/JinyuanSun/Codes_for_FoldX/blob/master/multiple_threads_foldx_positionscan.py) for details.  
If you have any questions, please contact 1650949260@qq.com


# PS:  
The FoldX5 provides a graphical user interface available as a python plugin for the YASARA molecular viewer.

# Reference
1. Javier Delgado, Leandro G Radusky, Damiano Cianferoni, Luis Serrano, FoldX 5.0: working with RNA, small molecules and a new graphical interface, Bioinformatics, Volume 35, Issue 20, 15 October 2019, Pages 4168–4169, https://doi.org/10.1093/bioinformatics/btz184
