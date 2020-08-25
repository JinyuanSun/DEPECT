#!/usr/bin/python3


import argparse
import os
import mutate_model_new

def get_parser():
    parser = argparse.ArgumentParser(description='active site single point mutation scan')
    #parser.parse_args()
    parser.add_argument("-l", "--ligand", help="ligand in PDBQT")
    parser.add_argument("-r", "--receptor", help="receptor in PDB")
    parser.add_argument("-bc", "--boxcfg", help="box config with box_size and box_center")
    parser.add_argument("-mf", "--mutation_file", 
                        help="single point mutations to do on the receptor structue, use with list mode")
    parser.add_argument("-pf", "--position_file", 
                        help="position to do saturated mutation on the receptor structue, use with scan mode")
    parser.add_argument("-m", "--mode", help="modes you want to use")
    parser.add_argument("-sf", "--scorefile", help="output file name",default="depect_enzde.sc")
    
    #args = parser.parse_args()
    return parser

#with pre-prepared ligand
#dock = single_dock.run_dock(name)

def vinalocalsearch(receptor,ligand,boxcfg):
    ADT_PATH = "/home/jsun/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/"
    preparereceptorcmd = "pythonsh "+ADT_PATH+"prepare_receptor4.py -r "+receptor+" -o "+receptor.replace(".pdb",".pdbqt")+" -A checkhydrogens"
    os.system(preparereceptorcmd)
    vinalocalcmd = "vina --receptor "+receptor.replace(".pdb",".pdbqt") +" --ligand "+ligand+" --local_only --config "+boxcfg+" --out "+ligand.replace(".pdbqt","@_"+receptor+"_out.pdbqt")
    AutoDockScoreCmd = "pythonsh "+ADT_PATH+"compute_AutoDock41_score.py -r "+receptor.replace(".pdb",".pdbqt")+" -l "+ligand.replace(".pdbqt","@_"+receptor+"_out.pdbqt")
    print(vinalocalcmd)
    localsearchoutput = os.popen(vinalocalcmd)
    os.system(AutoDockScoreCmd)
    for line in localsearchoutput.read().split("\n"):
        if line.startswith("Affinity:"):
            affinity = line.split(":")[1].split("(")[0].strip()
            return affinity

def outputfile(mutant,affinity,scorefilename):

    with  open(scorefilename,"a+") as scorefile:
        scorefile.write(mutant+","+str(affinity)+"\n")
        scorefile.close()
    
    

                
#_split("my_docking.pdb")
def _read_mutation_file(mut_file):
    #position|AA name|chain
    mut_list = []
    mut_file = open(mut_file)
    for line in mut_file:
        lst = line.strip().split("|")
        pac = [lst[0],lst[1],lst[2]]
        mut_list.append(pac)
    return mut_list

def _read_position_file(pos_file):
    #position|AA|chain
    pos_list = []
    pos_file = open(pos_file)
    for line in pos_file:
        lst = line.strip().split("|")
        pAc = [lst[0],lst[1],lst[2]]
        pos_list.append(pAc)
    return pos_list
        
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    receptor = args.receptor
    ligand = args.ligand
    mode = args.mode
    boxcfg = args.boxcfg
    scorefilename = args.scorefile
    ADT_PATH = "/home/jsun/MGLTools/MGLToolsPckgs/AutoDockTools/Utilities24/"
    with  open(scorefilename,"w+") as scorefile:
        scorefile.write("#mutant"+","+"vina_affinity"+"\n")
        scorefile.close()
    if mode == "list":
        mut_file = args.mutation_file
        mut_list = _read_mutation_file(mut_file)
        for pac in mut_list:
            print("Building mutation "+" ".join(pac))
            #mutate(pdb,respos,resname,chain):
            mutant = mutate_model_new.run_mutate(receptor,pac[0],pac[1],pac[2])
            outputfile(mutant,vinalocalsearch(mutant,ligand,boxcfg),scorefilename)
    if mode == "scan":
        pos_file = args.position_file
        pos_list = _read_position_file(pos_file)
        for pAc in pos_list:
            for X in ['CYS', 'ASP', 'SER', 'GLN', 'LYS','ILE', 'PRO', 'THR', 
                      'PHE', 'ASN','GLY', 'HIS', 'LEU', 'ARG', 'TRP','ALA', 'VAL', 'GLU', 'TYR', 'MET']:
                if pAc[1] != X:
                    mutant = mutate_model_new.run_mutate(receptor,pAc[0],X,pAc[2])
                    outputfile(mutant,vinalocalsearch(mutant,ligand,boxcfg),scorefilename)
