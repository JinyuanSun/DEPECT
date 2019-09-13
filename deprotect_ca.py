import numpy as np
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='prepare for BuildModel in FoldX')
parser.add_argument("-s", '--sequence', help="input a single protein sequence in fasta format!")
parser.add_argument("-i", '--identity cutoff', help='the identity cutoff when build the sub-database')
parser.add_argument("-e", '--evalue cutoff', help='the evalue cutoff when blast')


args = parser.parse_args()

seqname = args.sequence
evalue = args.evalue
idcutoff = int(args.identity)


dbname = "GT_60_reduced.fastadb"


def _build_hot_dic(dbname):
    hot_dic = {}
    hotdb = open(dbname)
    for line in hotdb:
        if line.startswith(">"):
            head = line.split(" ")[0][1:]
            hot_dic[head] = ""
        else:
            hot_dic[head] += line.replace("\n", "")
    return hot_dic


hot_dic = _build_hot_dic(dbname)


# 1. blastp against the therbase

def _blast(seqname, db, evalue, idcutoff):
    outfile = seqname + "_out.txt"
    cmd = "blastp -query " + seqname + " -db " + db + " -evalue " + evalue + " -outfmt 6 -out " + outfile
    # print(cmd)
    os.system(cmd)
    l = os.popen("wc -l " + outfile)
    ll = l.read()
    # print("blast finished, total "+ll+" thermostable sequences found")
    blastout = open(outfile)
    hmhits = []
    h = 0
    for line in blastout:
        w = line.replace("\n", "").split("\t")
        # print(w)
        identity = float(w[2])
        if identity > idcutoff:
            hmhits.append(w[1])
            h = h + 1
    # print("find "+str(h)+" homoglous sequences")
    subdb_name = seqname + ".prodb"
    os.system("cp " + seqname + " " + subdb_name)
    subdb = open(subdb_name, "a+")
    outline = ''
    for x in hmhits:  # extract seq from db for align
        # print(x)
        outline += ">" + x + "\n"
        outline += hot_dic.get(x) + "\n"
    print(outline, file = subdb)
    # os.system("cat "+seqname+" "+subdb_name+" > seq_hot.fasta")
    print("File " + subdb_name + " written!" + "\n" + "Alignment starts!")
    return subdb_name


# In[53]:


db = "GT_60_reduced.fastadb"
#evalue = "10"
#idcutoff = 30
#seqname = "petase.fasta"
subdb_name = _blast(seqname, db, evalue, idcutoff)

import os
from collections import Counter

seq_db_name = subdb_name


def _target(seqname):
    f = open(seqname)
    for l in f:
        if l.startswith(">"):
            target = l.replace("\n", "")
    return target


target = _target(seqname)


# target = ">sp|A0A0K8P6T7|PETH_IDESA Poly(ethylene terephthalate) hydrolase OS=Ideonella sakaiensis (strain 201-F6) OX=1547922 GN=ISF6_4831 PE=1 SV=1"
def align_analysis(seq_db_name, target):
    ali_filename = seq_db_name + ".ali"
    os.system("mafft" + "  --quiet  --localpair  --maxiterate 16  " + seq_db_name + " > " + ali_filename)
    ali_file = open(ali_filename)
    ali_dic = {}
    for line in ali_file:
        if line.startswith(">"):
            head = line.replace("\n", "")
            ali_dic[head] = ''
        else:
            seq = line.replace("\n", "")
            ali_dic[head] += seq
    # print(ali_dic)
    tseq = ali_dic.get(target)
    # print(tseq)
    n = 0
    p = 0
    pmlist = []
    for x in tseq:
        wt = x
        if wt == "-":
            n = n + 1
            # print(n)
        else:
            p = p + 1
            htlist = []
            for k in ali_dic:
                if k == target:
                    continue
                else:
                    hpm = ali_dic[k][n]
                    htlist.append(hpm)
            d = Counter(htlist)
            # print(d)
            l = len(htlist)
            # print(l)
            # if d.most_common()[0][1]/l < 0.3:
            # n = n +1
            # if d.most_common()[0][0] == "-":
            # n = n +1
            # if d.most_common()[0][0] == wt:
            # continue
            # else:
            pmlist.append([wt, p, d.most_common()])
            # print(wt,p,d.most_common()[0][0])
            n = n + 1
            # print(n)
        # print(n)
    return pmlist


pmlist = align_analysis(seq_db_name, target)


# print(pmlist)

def _get_num(subdb_name):
    n = 0
    f = open(subdb_name)
    for l in f:
        if l.startswith(">"):
            n = n + 1
    return n


l = _get_num(subdb_name)
# l = 43
cutoff = 0.5
for x in pmlist:
    wt = x[0]
    p = x[1]
    d = x[2]
    if d[0][0] == "-":
        try:
            mt = d[1][0]
            f = d[1][1]
            # print(wt,p,mt)
        except IndexError:
            continue
    else:
        mt = d[0][0]
        f = d[0][1]
    if f / l > cutoff:
        if mt == wt:
            continue
        else:
            print(wt, p, mt)
