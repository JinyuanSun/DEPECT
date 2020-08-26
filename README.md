# DEPECT
 Design and Engineering of Proteins and Enzymes by Computational Tools  
# Introduction
Enzymes are key to bioindustry and some proteins are very important like COVID-19 neutralize antibody. However, it is hard to get efficient or useful biomacromolecules. Here I present DEPECT package for more easily and faster in silico design of biomacromolecules.
# [enzde](https://github.com/JinyuanSun/DETECT/tree/master/enzde)
combine autodock suites and modeller to design enzymes.  
### To run enzde scripts, make sure you have AutoDock suites and modeller installed!  
 - Install  
 Edit the depect_enzde_config.py, change the ADT_PATH to where its installed.
 - Demo  
 `./demo.sh`

# Doc
## mutfile format
```
12|LYS|A  
|  |   |  
|  |   |___Chain ID  
|  |___Mutate to which AA (Here is to mutate to a LYS)   
|___Residue Sequence Number in Pdb File   
```

