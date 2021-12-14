#!/usr/bin/python3


import numpy as np
import os
import struct

"""
class Protein:
    def __init__(self, residue, atom):
        self.residue = residue
        self.atom = atom
    
    def get_atom_coords(self, residue, atom):
"""
# from time import time
# t0 = time()


def readpdb(pdbfilename):
    pdb_format = "6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s10s2s3s"
    protein_pdbdict = {}
    ligand_pdbdict = {}
    pdbfile = open(pdbfilename)
    for line in pdbfile:
        # try:
        # print(type(line))
        try:
            tmpline = struct.unpack(pdb_format, line.encode())
            # print(tmpline)
            if tmpline[0].decode().strip() == "ATOM":
                # record_name = line[0:6].replace(" ","")
                atom_num = tmpline[1].decode().strip()
                atom = tmpline[3].decode().strip()
                altLoc = tmpline[4].decode().strip()
                res = tmpline[5].decode().strip()
                chainid = tmpline[7].decode().strip()
                resseq = tmpline[8].decode().strip()
                icode = tmpline[10].decode().strip()
                x = float(tmpline[11].decode().strip())
                y = float(tmpline[12].decode().strip())
                z = float(tmpline[13].decode().strip())
                occ = tmpline[14].decode().strip()
                tfactor = tmpline[15].decode().strip()
                element = tmpline[17].decode().strip()
                charge = tmpline[18].decode().strip()
                try:
                    protein_pdbdict[chainid][resseq][atom] = np.array([x, y, z])
                except KeyError:
                    try:
                        protein_pdbdict[chainid][resseq] = {}
                        protein_pdbdict[chainid][resseq]["res"] = res
                        protein_pdbdict[chainid][resseq][atom] = np.array([x, y, z])
                    except KeyError:
                        protein_pdbdict[chainid] = {resseq: {atom: np.array([x, y, z])}}

            if tmpline[0].decode().strip() == "HETATM":
                # record_name = line[0:6].replace(" ","")
                atom_num = tmpline[1].decode().strip()
                atom = tmpline[3].decode().strip()
                altLoc = tmpline[4].decode().strip()
                res = tmpline[5].decode().strip()
                chainid = tmpline[7].decode().strip()
                resseq = tmpline[8].decode().strip()
                icode = tmpline[10].decode().strip()
                x = tmpline[11].decode().strip()
                y = tmpline[12].decode().strip()
                z = tmpline[13].decode().strip()
                occ = tmpline[14].decode().strip()
                tfactor = tmpline[15].decode().strip()
                element = tmpline[17].decode().strip()
                charge = tmpline[18].decode().strip()
                try:
                    ligand_pdbdict[chainid][atom_num][atom] = np.array([x, y, z])
                except KeyError:
                    try:
                        ligand_pdbdict[chainid][atom_num] = {}
                        ligand_pdbdict[chainid][atom_num]["res"] = res
                        ligand_pdbdict[chainid][atom_num][atom] = np.array([x, y, z])
                    except KeyError:
                        ligand_pdbdict[chainid] = {resseq: {atom: np.array([x, y, z])}}
        except struct.error:
            pdb_format = "6s5s1s4s1s3s1s1s4s1s3s8s8s8s6s6s10s2s1s"
            try:
                tmpline = struct.unpack(pdb_format, line.encode())
                if tmpline[0].decode().strip() == "ATOM":
                    # record_name = line[0:6].replace(" ","")
                    atom_num = tmpline[1].decode().strip()
                    atom = tmpline[3].decode().strip()
                    altLoc = tmpline[4].decode().strip()
                    res = tmpline[5].decode().strip()
                    chainid = tmpline[7].decode().strip()
                    resseq = tmpline[8].decode().strip()
                    icode = tmpline[10].decode().strip()
                    x = float(tmpline[11].decode().strip())
                    y = float(tmpline[12].decode().strip())
                    z = float(tmpline[13].decode().strip())
                    occ = tmpline[14].decode().strip()
                    tfactor = tmpline[15].decode().strip()
                    element = tmpline[17].decode().strip()
                    charge = tmpline[18].decode().strip()
                    try:
                        protein_pdbdict[chainid][resseq][atom] = np.array([x, y, z])
                    except KeyError:
                        try:
                            protein_pdbdict[chainid][resseq] = {}
                            protein_pdbdict[chainid][resseq]["res"] = res
                            protein_pdbdict[chainid][resseq][atom] = np.array([x, y, z])
                        except KeyError:
                            protein_pdbdict[chainid] = {
                                resseq: {atom: np.array([x, y, z])}
                            }
            except struct.error:
                # print(line)
                continue
    return {"protein": protein_pdbdict, "ligand": ligand_pdbdict}


def readpdbqt(pdbqtfilename):
    pdbqtfile = open(pdbqtfilename)
    pdbqt_dic = {}
    for line in pdbqtfile:
        try:
            tmp = line.split()
            type = tmp[0]
            atom_num = tmp[1]
            atom = tmp[2]
            lig = tmp[3]
            coordx = float(tmp[5])
            coordy = float(tmp[6])
            coordz = float(tmp[7])
            if type == "HETATM":
                pdbqt_dic[atom_num] = np.array([coordx, coordy, coordz])
        except IndexError:
            continue
        except ValueError:
            continue
    return pdbqt_dic


# for i in range(100):
# a = readpdb("example/4hs9.pdb")
# print("pdb2dict parse a pdb",time()-t0)

# protein_pdbdict = readpdb("4hs9_p_A_ARG_12.pdb")["protein"]
# ligand_pdbdict = readpdbqt("example/RPBE_ref_out.pdbqt")
# ["A"]["79"]["OG"])
# print(protein_pdbdict)
# print(protein_pdbdict["protein"]["A"]["79"]["OG"])
# p_A_79_OG = protein_pdbdict["A"]["79"]["OG"]
# l_13_O = ligand_pdbdict["13"]
# from scipy.spatial import distance
# a = (1, 2, 3)
# b = (4, 5, 6)
# dst = distance.euclidean(l_13_O, p_A_79_OG)
# print(type(l_13_O[0]))
###############################################
# change the ligand name


def caldistance(ligand, depectcafile):
    with open("depect_enzde_distance.sc", "w+") as newoutfile:
        newoutfile.write("#receptor,vina_affinity,distance1\n")
        newoutfile.close()
    for line in depectcafile:
        if line.startswith("#"):
            continue
        else:
            tmp = line.split(",")
            receptor = proteinfile = tmp[0]
            ligandfile = ligand.replace(".pdbqt", "@_" + receptor + "_out.pdbqt")
            # ligandfile = ligand.replace(".pdbqt",receptor+"_out.pdbqt")
            protein_pdbdict = readpdb(proteinfile)
            ligand_pdbdict = readpdbqt(ligandfile)
            # ["A"]["79"]["OG"])
            # print(protein_pdbdict["protein"]["A"]["79"]["OG"])
            p_A_79_OG = protein_pdbdict["protein"]["A"]["79"]["OG"]
            l_13_O = ligand_pdbdict["13"]
        # print(proteinfile,ligandfile)
        distance = round(np.linalg.norm(p_A_79_OG - l_13_O), 3)
        with open("depect_enzde_distance.sc", "a+") as newoutfile:
            newoutfile.write(line.replace("\n", "," + str(distance) + "\n"))
            newoutfile.close()


if __name__ == "__main__":
    ligand = "RPBE_ref.pdbqt"
    depectcafile = open("depect_enzde.sc")
    caldistance(ligand, depectcafile)
