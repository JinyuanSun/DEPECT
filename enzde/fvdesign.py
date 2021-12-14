#!/usr/bin/env python

#
import os
import time
import numpy as np
import pandas as pd
import argparse
from joblib import Parallel, delayed


class FoldX:
    def __init__(self, exe=""):
        if exe:
            self.exe = exe
        else:
            self.exe = os.popen("which foldx").read().replace("\n", "")

        # AlaScan
        # AnalyseComplex
        # BuildModel
        # CrystalWaters
        # Dihedrals
        # DNAContact
        # DNAScan
        # RNAScan
        # MetalBinding
        # Optimize
        # PDBFile
        # PositionScan
        # PrintNetworks
        # Pssm
        # QualityAssessment
        # ReconstructSideChains
        # RepairPDB
        # Rmsd
        # SequenceDetail
        # SequenceOnly
        # Stability

    def repair(self, input_file):
        cmd_list = [self.exe, "--command=RepairPDB", "--pdb=%s" % input_file]

    def build_model(self, input_file, mutation, numberOfRuns=1):
        # mutation = "A_126_M"
        start_time = time.time()
        cmd_list = [
            self.exe,
            "--command=BuildModel",
            "--pdb=%s" % input_file,
            "--mutant-file=individual_list.txt",
            "--numberOfRuns=%s 1>/dev/null" % str(numberOfRuns),
        ]
        with open("individual_list.txt", "w+") as indifile:
            indifile.write(mutation + ";\n")
            indifile.close()
        # os.popen(" ".join(cmd_list))
        if test_on_mac:
            print(" ".join(cmd_list))
        else:
            os.system(" ".join(cmd_list))
        end_time = time.time()
        return round(end_time - start_time, 3)


class Autodock:
    def __init__(self, pythonsh_path, ADTU_path, box_cache, ligand_cache):
        # /home/jsun/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24
        self.ADTU_path = ADTU_path
        self.pythonsh_path = pythonsh_path
        self.box_cache = box_cache
        self.ligand_cache = ligand_cache

    def run_local_dock(self, receptor, i):
        os.system(
            "%s %s/prepare_receptor4.py -r %s -A checkhydrogens"
            % (self.pythonsh_path, self.ADTU_path, receptor)
        )
        # os.popen("pythonsh %s/prepare_ligand.py -l %s " %(self.ADTU_path, ligand))
        with open("box.cfg", "w+") as box:
            box.write(self.box_cache)
            box.close()
        with open("ligand.pdbqt", "w+") as ligand:
            ligand.write(self.ligand_cache)
            ligand.close()
        local_out = os.popen(
            "vina --config box.cfg --receptor=%s --ligand=ligand.pdbqt --local_only --cpu=1 --out=local_docked_%s.pdbqt"
            % (receptor + "qt", str(i))
        )
        for line in local_out.read().split("\n"):
            if line.startswith("Affinity:"):
                affinity = line.split(":")[1].split("(")[0].strip()
                return affinity


def _3_2_1(x):
    d = {
        "CYS": "C",
        "ASP": "D",
        "SER": "S",
        "GLN": "Q",
        "LYS": "K",
        "ILE": "I",
        "PRO": "P",
        "THR": "T",
        "PHE": "F",
        "ASN": "N",
        "GLY": "G",
        "HIS": "H",
        "LEU": "L",
        "ARG": "R",
        "TRP": "W",
        "ALA": "A",
        "VAL": "V",
        "GLU": "E",
        "TYR": "Y",
        "MET": "M",
    }
    assert x in d, "%s is not in 20 canonical amino acids!" % (x)
    return d[x]


def __remove_chain_from_individual_list_mutation(mutation):
    # EB243Q -> [[E,B,243,Q]]
    # EA243Q,EB243Q -> [[E,A,243,Q],[E,B,243,Q]]
    mutation_list = mutation.split(",")
    converted = []
    for mut in mutation_list:
        wild_type = mut[0]
        mut_type = mut[-1]
        chain = mut[1]
        num = mut[2:-1]
        converted.append([wild_type, chain, num, mut_type])
    return converted


class Pipeline:
    def __init__(
        self,
        position_list_filename,
        enzyme_wildtype_structure_filename,
        substrate_filename,
        chain_list,
        numberOfRuns,
        box_cfg,
        ADTU_path="",
        pythonsh_path="",
        verbose=False,
    ):
        self.text = 0
        self.position_list_filename = position_list_filename
        self.chain_dict = {}
        self.root_dir = os.getcwd()
        self.enzyme = "/".join([self.root_dir, enzyme_wildtype_structure_filename])
        self.substrate = "/".join([self.root_dir, substrate_filename])
        self.lig = open(self.substrate, "r").read()
        self.verbose = verbose
        self.cachefile = open(enzyme_wildtype_structure_filename, "r").read()
        # self.docking_score = {}
        self.box_cache = open(box_cfg, "r").read()
        self.numberOfRuns = numberOfRuns
        self.chain_list = chain_list
        # self.foldx_energy = {}

        if ADTU_path:
            self.ADTU_path = ADTU_path
        else:
            self.ADTU_path = (
                os.popen("locate MGLToolsPckgs/AutoDockTools/Utilities24")
                .read()
                .split("\n")[0]
            )

        if pythonsh_path:
            self.pythonsh_path = pythonsh_path
        else:
            self.pythonsh_path = os.popen("which pythonsh").read().replace("\n", "")

    def read_in_positions(self):
        # mutation list generated from "pymol selection sele, output=S"
        with open(self.position_list_filename, "r") as pos_file:
            for line in pos_file:
                if line[0] != "#":
                    protein_name, chain, wild_type_3l, res_num = line.replace(
                        "\n", ""
                    ).split("/")
                    # self.position_list.append([chain, wild_type, res_num])
                    if chain in self.chain_dict:
                        self.chain_dict[chain].append((wild_type_3l, res_num))
                    else:
                        self.chain_dict[chain] = [(wild_type_3l, res_num)]
            pos_file.close()

    def generate_all_mutation(self):
        # Automate determine how many chains to use depending on input list file
        docking_score = {}
        foldx_energy = {}
        if self.chain_list:
            chains = self.chain_list.split(",")
        else:
            chains = list(self.chain_dict.keys())
        poses = []

        for chain in self.chain_dict:
            poses += self.chain_dict[chain]
        all_mutation = []

        for pos in poses:
            pos_mutation_info = []
            wild_type_3l, res_num = pos
            wild_type_1l = _3_2_1(wild_type_3l)
            for mut_type in "QWERTYIPASDFGHKLCVNM":
                if mut_type != wild_type_1l:
                    a_mutation = []
                    for chain in chains:
                        a_mutation.append(
                            "".join([wild_type_1l, chain, res_num, mut_type])
                        )
                    pos_mutation_info.append(",".join(a_mutation))
                    docking_score[",".join(a_mutation)] = np.zeros(
                        self.numberOfRuns
                    )
                    foldx_energy[",".join(a_mutation)] = {}
            all_mutation += pos_mutation_info

        return all_mutation, docking_score, foldx_energy

    def convert_mut_file_to_mutations(self):
        pass

    def calScore(self, mutation, foldx_energy):
        # fxout_name = jobID + "/Dif_" + pdbfile.replace(".pdb", ".fxout")
        fxout_name = os.popen("ls Dif*.fxout").read().replace("\n", "")
        # print(fxout_name)
        df = pd.read_table(fxout_name, sep="\t", skiprows=8)
        score = round(df["total energy"].mean(), 4)
        sd = round(df["total energy"].std(), 4)
        foldx_energy[mutation] = {"score": score, "SD": sd}
        with open("foldx_energy.txt", 'w+') as output:
            output.write("\t".join([mutation, str(score), str(sd)]))
            output.close()
        # return ["_".join([wild, str(resNum), mutation]), score, sd]

    def build_model(self, varlist):
        mutation, docking_score, foldx_energy = varlist
        os.mkdir(mutation)
        os.chdir(mutation)
        os.mkdir("utils/build_model")
        os.chdir("utils/build_model")
        with open("WT_protein.pdb", "w+") as enzymefile:
            enzymefile.write(self.cachefile)
            enzymefile.close()
        runtime = FoldX().build_model("WT_protein.pdb", mutation, self.numberOfRuns)
        self.calScore(mutation, foldx_energy)

        if self.substrate:
            i = 0
            aff_arr = np.zeros(self.numberOfRuns)
            while i in range(self.numberOfRuns):
                affinity = Autodock(
                    pythonsh_path=self.pythonsh_path,
                    ADTU_path=self.ADTU_path,
                    box_cache=self.box_cache,
                    ligand_cache=self.lig,
                ).run_local_dock("WT_protein_1_%s.pdb" % str(i), i)
                # docking_score[mutation][i] = affinity
                aff_arr[i] = affinity
                i += 1
            with open("docking_energy.txt", 'w+') as output:
                output.write("\t".join([mutation, str(round(aff_arr.mean(), 4)), str(round(aff_arr.std(), 4))]))
        os.chdir("../")
        if self.verbose:
            print("[DEBUG]: ")

    def build_scan(self, mutations, docking_score, foldx_energy, threads):

        Parallel(n_jobs=threads)(
            delayed(self.build_model)([mutation, docking_score, foldx_energy]) for mutation in mutations
        )

        return 0


def get_args():
    parser = argparse.ArgumentParser(
        description="Run FoldX, AutoDock Vina for substrate binding pocket redesign."
    )
    parser.add_argument("enzyme", help="Input enzyme PDB")
    parser.add_argument(
        "chain_list", help='list of chains, capital letters seperated with comma: "A,B"'
    )
    parser.add_argument("substrate", help="Input substrate PDBQT")
    parser.add_argument("box_cfg", help="box information for local docking")
    parser.add_argument("position_list", help="position list file output from pymol")
    parser.add_argument("-T", "--threads", help="Number of threads to use", default=16)
    parser.add_argument(
        "-N",
        "--num_of_runs",
        help="Number of model foldx to build for each mutation",
        default=5,
    )
    parser.add_argument("-mp", "--mgltools_path", help="path to mgltools")

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = get_args()
    test_on_mac = False
    position_list_filename = args.position_list
    enzyme_wildtype_structure_filename = args.enzyme
    substrate_filename = args.substrate
    box_cfg = args.box_cfg
    chain_list = args.chain_list
    numberOfRuns = int(args.num_of_runs)
    threads = int(args.threads)
    mgltools_path = args.mgltools_path
    scan_pipeline = Pipeline(
        position_list_filename,
        enzyme_wildtype_structure_filename,
        substrate_filename,
        chain_list,
        numberOfRuns,
        box_cfg,
        mgltools_path,
    )

    def scan_and_dock():
        scan_pipeline.read_in_positions()
        all_mutations, docking_score, foldx_energy = scan_pipeline.generate_all_mutation()
        scan_pipeline.build_scan(all_mutations[:6], docking_score, foldx_energy, threads)

    scan_and_dock()
