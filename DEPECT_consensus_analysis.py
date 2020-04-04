#######################################################
#                      DEPECT                         #
#                                                     #
#    DEsign Protein and Enzyme by Computation Tools   #
#                                                     #
#                Author: Sun Jinyuan                  #
#                                                     #
#######################################################

import argparse
import os
import math
import time


def get_parser():
    parser = argparse.ArgumentParser(description='consensus analysis based on relative entropy score')
    #parser.parse_args()
    parser.add_argument("-s", "--sequence", help="target sequence which must be provided")
    parser.add_argument("-db", "--database", help="database name for blastp")
    parser.add_argument("-ic", '--identity_cutoff',default=30,
                        help='the identity cutoff when build the sub-database')
    parser.add_argument("-e", '--evalue', default=1e-5,
                        help='the evalue cutoff when blast')
    parser.add_argument("-ds", "--dataset", help="dataset to do consensus analysis")
    parser.add_argument("-sc", '--score_cutoff', default=2.25,
                        help='the score cutoff for selection, default is 2.25')
    parser.add_argument("-m", '--mode', help='blast mode require input -ic, -db and -e, analysis mode require -ds')
    #args = parser.parse_args()
    return parser

def build_sub_db(seq,db_name,identity_cutoff,evalue):
    blastout = os.popen("blastp -query "+seq+" -db "+str(db_name)+" -evalue "+str(evalue)+" -outfmt 6")
        #split blastout into list
    try:
        blastout = blastout.read()[:-2].split("\n")
        with open("hits_list.txt","w+") as hitlist:
            for line in blastout:
                identity = float(line.split("\t")[2])
                hit = line.split("\t")[1]
                if identity > int(identity_cutoff):
                    hitlist.write(hit+"\n")
            hitlist.close()
        #hitlist = open("hits_list.txt","w+")
        #hitlist = open("hits_list.txt","a+")
            
        os.system("blastdbcmd -db "+db_name+" -entry_batch hits_list.txt -out target_refdb.fasta")
        #refdb = os.popen("blastdbcmd -db "+db_name+" -entry_batch hits_list.txt")
        #refdb = refdb.read()
        #refdb_file = open("target_refdb.fasta","w+")
        #refdb_file = open("target_refdb.fasta","a+")
        #print(refdb,file=refdb_file)
        os.system("cat "+seq+" target_refdb.fasta > "+seq+".db")
        os.system("muscle -in "+seq+".db -out "+seq+".aln -maxiters 99 -quiet")
        return seq+".aln"
    except IndexError:
        print("Can't find any hit, please change evalue or identity cutoff")
        return "No match"

def cal_score(aa_list):
    length = len(aa_list)
    aa_count_dict = {}
    aa_score_dict = {}
    D = 0
    for aa in aa_list:
        if aa in aa_count_dict:
            aa_count_dict[aa] += 1
        else:
            aa_count_dict[aa] = 1
    AA_composition={"A":  8.25, "R":  5.53, "N":  4.06, "D":  5.45, "C":  1.37, 
                "Q":  3.93, "E":  6.75, "G":  7.07, "H":  2.27, "I":  5.96,
                "L":  9.66, "K":  5.84, "M":  2.42, "F":  3.86, "P":  4.70,
                "S":  6.56, "T":  5.34, "W":  1.08, "Y":  2.92, "V":  6.87}
    for aa in aa_count_dict:
        if aa != "-":
            aa_p = aa_count_dict[aa]/length
            aa_q = AA_composition[aa]/100
            D = D + aa_p*math.log(aa_p/aa_q)
    for aa in aa_count_dict:
        if aa != "-":
            aa_p = aa_count_dict[aa]/length
            aa_score = aa_p*D
            aa_score_dict[aa] = aa_score
    return aa_score_dict 
    
def build_matrix(aln,seq):
    map_dict = {}
    matrix_dict = {}
    score_dict = {}
    head_alnseq_dict = {}
    seq_file = open(seq)
    raw_target = ""
    #read target name
    for line in seq_file:
        if line.startswith(">"):
            target_name = line.strip()
        else:
            raw_target = raw_target + line.strip()
    aln_file = open(aln)
    #build head_alnseq_dict
    for line in aln_file:
        if line.startswith(">"):
            head = line.strip()
            head_alnseq_dict[head] = ""
        else:
            head_alnseq_dict[head] += line.strip()
    #build map_dict to set the position relation between the aligned target and raw sequence
    aln_target = head_alnseq_dict[target_name]
    raw_pos = 0
    aln_pos = 0
    for aa in aln_target:
        if aa == "-":#a gap
            raw_pos = raw_pos
            aln_pos = aln_pos + 1
            map_dict[aln_pos] = raw_pos
        if aa != "-":#an amino acid
            raw_pos = raw_pos + 1
            aln_pos = aln_pos + 1
            map_dict[aln_pos] = raw_pos
    #build matrix_dict
    for key in head_alnseq_dict:
        alnseq = head_alnseq_dict[key]
        n = 0
        for aa in alnseq:
            n = n + 1
            try:
                matrix_dict[n].append(aa)
            except KeyError:
                matrix_dict[n] = [aa]
    for position in matrix_dict:
        for aa_list in matrix_dict[position]:
            aa_score_dict = cal_score(aa_list)
            score_dict[position] = aa_score_dict
    return score_dict,map_dict,raw_target

def select(score_dict,map_dict,score_cutoff,raw_target):
    out_line = 'position\tmutation\tscore'
    for aln_pos in score_dict:
        aa_score_dict=score_dict[aln_pos]
        for aa in aa_score_dict:
            if aa != "-":
                score = aa_score_dict[aa]
                if score > float(score_cutoff):
                    raw_pos = map_dict[aln_pos]
                    raw_aa = raw_target[raw_pos-1]
                    if aa != raw_aa:
                        out_line = out_line + str(raw_pos) +"\t"+ aa +"\t"+ str(score)[0:4] + "\n"
                    
    return out_line



if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    seq = args.sequence
    db_name = args.database
    identity_cutoff = args.identity_cutoff
    evalue = args.evalue
    score_cutoff = args.score_cutoff
    mode = args.mode
    if mode == "blast":
        aln = build_sub_db(seq,db_name,identity_cutoff,evalue)
        score_dict,map_dict,raw_target = build_matrix(aln,seq)
        out_line = select(score_dict,map_dict,score_cutoff,raw_target)
        print(out_line)
    if mode == "analysis":
        dataset = args.dataset
        os.system("cat "+seq+" "+dataset+" > "+seq+".db")
        os.system("muscle -in "+seq+".db -out "+seq+".aln -maxiters 99 -quiet")
        aln = seq+".aln"
        score_dict,map_dict,raw_target = build_matrix(aln,seq)
        out_line = select(score_dict,map_dict,score_cutoff,raw_target)
        print(out_line)
