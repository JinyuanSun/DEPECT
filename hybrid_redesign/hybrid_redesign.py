#!/usr/bin/python3
import os
import math
import argparse
import time

import depect_ca_local as dcl
#import rosetta_cartddg_filter
import ca2mutfile

#cartesian_ddg.linuxgccrelease -s 4eb0_clean.pdb @cart_ddg_flag -ddg:mut_file 4EB01.mutfile
Cartesian_ddG_exe = "cartesian_ddg.linuxgccrelease"
Cartesian_ddG_opt = "-ddg:iterations 1 -ddg::cartesian -ddg::dump_pdbs true -ddg:bbnbrs 1 -fa_max_dis 9.0 -score:weights ref2015_cart"


def get_parser():
    parser = argparse.ArgumentParser(description='consensus analysis based on relative entropy score and fruther filter by ddg')
    #parser.parse_args()
    parser.add_argument("-s", "--structure", help="target structure which must be provided")
    parser.add_argument("-cid", "--chain", help="chain id for sequence extraction, default is A",default = "A")
    parser.add_argument("-db", "--database", help="database name for blastp")
    parser.add_argument("-ic", '--identity_cutoff',default=30,
                        help='the identity cutoff when build the sub-database')
    parser.add_argument("-e", '--evalue', default=1e-5,
                        help='the evalue cutoff when blast')
    parser.add_argument("-ds", "--dataset", help="dataset to do consensus analysis")
    parser.add_argument("-sc", '--score_cutoff', default=2,
                        help='the score cutoff for selection, default is 2')
    parser.add_argument("-nt", '--num_threads', default=8,
                        help='number of threads used to run blastp and ddg calculation software, default is 8')
    parser.add_argument("-m", '--mode', help='blast mode require input -ic, -db and -e, analysis mode require -ds')
    parser.add_argument("-en", '--engine', default="rcd",
                        help='software used for ddg calculation, currently supporting FoldX (f) and Rosetta_Cartesian_ddG (rcd). Default is Rosetta_Cartesian_ddG')
    parser.add_argument("-dc", '--ddg_cutoff', help='cutoff of ddg for output')
    #args = parser.parse_args()
    return parser

def pdb2seq(structure,chain):
    pdbfile = open(structure)
    seq = ''
    #pdb_seq_out_dict = {}
    longer_names = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                  'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                  'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                  'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                  'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
    #pdb = open(pdbname)
    resseqlst = []
    for line in pdbfile:
        if line.startswith("ATOM") and "CA" in line.split():
            amino_acid = line[17:20]
            resseq = int(line[22:26].replace(" ",""))
            if chain == line[21:22] and resseq not in resseqlst:
                resseqlst.append(resseq)
                seq = seq + longer_names[amino_acid]
        #if line[0:3] == "TER":
            #seq = seq + "\n"
    #pdb_seq_out_dict[chain] = seq[:-1]
    seqfilename = structure.replace("pdb","fasta")
    with open(seqfilename,"w+") as seqfile:
        seqfile.write(">seq\n"+seq+"\n")
        seqfile.close()
    return seqfilename

def runconsensusanalysis(seq,db_name,identity_cutoff,evalue,nt,score_cutoff,mode):
    
    #parser = dcl.get_parser()
    #args = parser.parse_args()
    #struct = args.structure
    #chain = args.chain
    #seq = pdb2seq(struct,chain)
    ##seq = args.structure.replace("pdb","fasta")
    #db_name = args.database
    #identity_cutoff = args.identity_cutoff
    #evalue = args.evalue
    #score_cutoff = args.score_cutoff
    #mode = args.mode
    #nt = str(args.num_threads)
    if mode == "blast":
        t0 =time.time()
        aln = dcl.build_sub_db(seq,db_name,identity_cutoff,evalue,nt)
        print("Blast search finished in",str(time.time()-t0)[0:3]+"s")
        score_dict,map_dict,raw_target = dcl.build_matrix(aln,seq)
        outlist = dcl.select(score_dict,map_dict,score_cutoff,raw_target)
        curtime = time.asctime( time.localtime(time.time()))
        depectcafile = seq+".ca"
        with open(depectcafile,"w+") as outf:
            outf.write("#Back-to-consensus mutations suggested by DEPECT consensus analysis.\n#"+curtime+"\n#wildtype\tnumber\tback-to-consensus\tscore"+"\n")
            outf.close()
        with open(depectcafile,"a+") as outf:
            for mutations in outlist:
                outf.write(mutations+"\n")
            outf.close()
        #print(out_line)
    if mode == "analysis":
        dataset = args.dataset
        os.system("cat "+seq+" "+dataset+" > "+seq+".db")
        os.system("muscle -in "+seq+".db -out "+seq+".aln -maxiters 99")
        aln = seq+".aln"
        score_dict,map_dict,raw_target = dcl.build_matrix(aln,seq)
        outlist = dcl.select(score_dict,map_dict,score_cutoff,raw_target)
        curtime = time.asctime( time.localtime(time.time()))
        depectcafile = seq+".ca"
        with open(depectcafile,"w+") as outf:
            outf.write("#Back-to-consensus mutations suggested by DEPECT consensus analysis.\n#"+curtime+"\n#wildtype\tnumber\tback-to-consensus\tscore"+"\n")
            outf.close()
        with open(depectcafile,"a+") as outf:
            for mutations in outlist:
                outf.write(mutations+"\n")
            outf.close()
    return depectcafile

def generaterosettamutfile(depectcafile,nt):
    #parser = dcl.get_parser()
    #args = parser.parse_args()
    nt = int(nt)
    mutfilelist = ca2mutfile.ca2mutfilelist(depectcafile)
    #print(nt)
    outfilename = prefix = depectcafile.split(".")[0]
    if nt == 1:
        ca2mutfile.mutfilelist2rosettamutfile(mutfilelist,outfilename)
    if nt > 1:
        mutfilenamelist = ca2mutfile.chunk(mutfilelist,nt,prefix)
    return mutfilenamelist

def runrosetta(mutfilenamelist,struct):
    #parser = dcl.get_parser()
    #args = parser.parse_args()
    #struct = args.structure
    for mutfile in mutfilenamelist:
        os.system("nohup "+Cartesian_ddG_exe+" -s "+struct+" "+Cartesian_ddG_opt+" -ddg:mut_file "+mutfile+" &")


if __name__ == '__main__':
    print(
    """
        #########################################################
        #                         DEPECT                        #
        #           thermostability  hybrid redesign            #
        #    Design and Engineering of Proteins and Enzymes     #
        #                by Computational Tools                 #
        #                                                       #
        #                  Author: Sun Jinyuan                  #
        #             E-mail: jinyuansun98@gmail.com            #
        #########################################################

        """
    )
    parser = get_parser()
    parser = get_parser()
    args = parser.parse_args()
    struct = args.structure
    chain = args.chain
    seq = pdb2seq(struct,chain)
    #seq = args.structure.replace("pdb","fasta")
    db_name = args.database
    identity_cutoff = args.identity_cutoff
    evalue = args.evalue
    score_cutoff = args.score_cutoff
    mode = args.mode
    nt = str(args.num_threads)
    depectcafile = runconsensusanalysis(seq,db_name,identity_cutoff,evalue,nt,score_cutoff,mode)
    mutfilenamelist = generaterosettamutfile(depectcafile,nt)
    runrosetta(mutfilenamelist,struct)
