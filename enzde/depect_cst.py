#!/usr/bin/python3

#
import math
import numpy as np
import os
import struct


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


def readcstfile(cstfilename):
    cstfile = open(cstfilename)
    cstparadict = {}
    # liganddict = readpdbqt(ligandfilename)
    # proteindict = readpdb(proteinfilename)
    for line in cstfile:
        atomlst = []
        tmp = line.replace("\n", "").split(":")
        head = tmp[0]
        atoms = tmp[1].split(",")
        for atom in atoms:
            atmp = atom.split("|")
            if atmp[0] == "protein":
                # coordlst.append(proteindict["protein"][atmp[1]][atmp[2]][atmp[3]])
                atomlst.append(["protein", atmp[1], atmp[2], atmp[3]])
            if atmp[0] == "ligand":
                atomlst.append(["ligand", atmp[1]])
        cstparadict[head] = atomlst
    return cstparadict


def cal_dihedral(p):
    """Praxeolitic formula
    1 sqrt, 1 cross product"""
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)
    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def getcoords(cstparadict, proteindict, liganddict):
    cstcoorddict = {}
    for head in cstparadict:
        atomlst = cstparadict[head]
        coordlst = []
        for atom in atomlst:
            if atom[0] == "protein":
                coordlst.append(proteindict["protein"][atom[1]][atom[2]][atom[3]])
            if atom[0] == "ligand":
                coordlst.append(liganddict[atom[1]])
        cstcoorddict[head] = coordlst
    return cstcoorddict


def calgeovalue(cstcoorddict, proteindict, liganddict):
    geovaldict = {}
    for key in cstcoorddict:
        coordlst = cstcoorddict[key]
        # print(coordlst)
        if len(coordlst) == 2:
            geovaldict[key] = str(round(np.linalg.norm(coordlst[0] - coordlst[1]), 3))
        if len(coordlst) == 3:
            geovaldict[key] = str(
                round(
                    math.acos(
                        getcos(coordlst[0] - coordlst[1], coordlst[1] - coordlst[2])
                    )
                    * 180
                    / math.pi,
                    3,
                )
            )
        if len(coordlst) == 4:
            geovaldict[key] = str(round(cal_dihedral(coordlst), 3))
    return geovaldict


def run_geo_cst(depectfilename, cstfilename, ligand):
    depectfile = open(depectfilename)
    cstparadict = readcstfile(cstfilename)
    cstfile = open(cstfilename)
    # liganddict = readpdbqt(ligandfilename)
    head_line = ""
    for line in cstfile:
        head_line = head_line + "," + line.split(":")[0]
    for line in depectfile:
        if line.startswith("#"):
            with open("depect_enzde_geo.sc", "w+") as newoutfile:
                newoutfile.write("#receptor,vina_affinity" + head_line + "\n")
                newoutfile.close()
        else:
            tmp = line.split(",")
            receptor = proteinfilename = tmp[0]
            ligandfilename = ligand.replace(".pdbqt", "@_" + receptor + "_out.pdbqt")
            # print(ligandfilename)
            # ligandfile = ligand.replace(".pdbqt",receptor+"_out.pdbqt")
            proteindict = readpdb(proteinfilename)
            liganddict = readpdbqt(ligandfilename)
            cstcoorddict = getcoords(cstparadict, proteindict, liganddict)
            # print(cstcoorddict)
            geovaldict = calgeovalue(cstcoorddict, proteindict, liganddict)
            # print(geovaldict)

            with open("depect_enzde_geo.sc", "a+") as newoutfile:
                newoutfile.write(
                    line.replace("\n", "," + ",".join(list(geovaldict.values())) + "\n")
                )
                newoutfile.close()
    # cstparadict = readcstfile(cstfilename,ligandfilename,proteinfilename)
    # geovaldict = calgeovalue(cstparadict)


if __name__ == "__main__":
    import sys

    try:
        depectfilename = sys.argv[1]
        cstfilename = sys.argv[2]
        ligand = sys.argv[3]
        run_geo_cst(depectfilename, cstfilename, ligand)
    except IndexError:
        print("Usage: python3 depect_cst.py depect_enzde.sc zlj.cst RPBE_ref.pdbqt")
        # try:
        # depectfilename = sys.argv[1]
        # cstfilename = sys.argv[2]
        # ligand = sys.argv[3]
        """
        try:
            open(depectfilename)
        except FileNotFoundError:
            print("Trying to run demo!")
            depectfilename = "depect_enzde.sc"
            cstfilename = "zlj.cst"
            ligand = "RPBE_ref.pdbqt"
            run_geo_cst(depectfilename,cstfilename,ligand)
            """
