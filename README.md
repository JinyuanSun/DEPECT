# DEPECT
 Design and Engineering of Proteins and Enzymes by Computational Tools  
# Introduction
Enzymes are key to bioindustry and some proteins are very important like COVID-19 neutralize antibody. However, it is hard to get efficient or useful biomacromolecules. Here I present DEPECT package for more easily and faster in silico design of biomacromolecules.
# [hybrid_redesign](https://github.com/JinyuanSun/DEPECT/tree/master/hybrid_redesign)
Click the [link](https://github.com/JinyuanSun/DEPECT/releases/download/v1.0/hybrid_redesign) to download the linux executable binary file of hybrid_redesign. It has been tested under Ubuntu 20.04.  
### To run hybrid_redesign scripts, make sure you have Rosetta, muscle and ncbi-blast+ suite installed!
#### usage:  
```
./hybrid_redesign -s 4eb0.pdb -cid A -db /ndata/databases/blastdb/nr -ic 30 -e 1e-20 -sc 2 -nt 16 -m blast
```
After the calculation done(You can monitor it via `top`), run `hrd.sh` to generate combined `hybrid_design.sc` file.

# [enzde](https://github.com/JinyuanSun/DETECT/tree/master/enzde)
combine autodock suites and modeller to design enzymes.  
### To run enzde scripts, make sure you have AutoDock suites and modeller installed!  
 - Install  
 Edit the depect_enzde_config.py, change the ADT_PATH to where its installed.
 - Demo  
 `./demo.sh`

# Doc
#### mutfile format
```
12|LYS|A  
|  |   |  
|  |   |___Chain ID  
|  |___Mutate to which AA (Here is to mutate to a LYS)   
|___Residue Sequence Number in Pdb File   
```
#### cstfile format
```
CST_NAME:atom1,atom2
If the atom in the protein:
protein|chain id|number of residue|atom name
If the atom in the ligand:
ligand|atom number
```
in the zlj.cst used in demo:  
```
distace1:protein|A|79|OG,ligand|13  
```
If the constraint is an angle, use 3 atoms, for a dihedral use 4 atoms.  
