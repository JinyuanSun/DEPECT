#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#for thermostable protein engineering
#learning knowlewdge from the nature
#how to evolve
#1. blastp against the therbase
#2. align with similiar sequences
#3. analysis


# In[1]:


import os
import numpy as np
import pandas as pd


# In[ ]:


#1. blastp against the therbase
db = ""
evalue = "1e-5"
cutoff = 30
seqname = ""
def _blast(seqname,db,evalue,cutoff):
    outfile = seqname+"_out.txt"
    cmd = "blastp -query "+seqname+" -db "+db+" -evalue "+evalue+" -outfmt 6 -out "+outfile
    os.system(cmd)
    l = os.popen("wc -l "+outfile)
    ll = l.read()
    print("blast finished, total "+ll+" thermostable sequences found")
    blastout = read(outfile)
    hmhits = []
    h = 0
    for line in blastout:
        w = line.replace("\n","").split("\t")
        identity = float(w[2])
        if identity > cutoff:
            hmhits += w[1]
            h = h + 1
    print("find "+h+" homoglous sequences")
    subdb_name = seqname+".prodb"
    subdb = open(subdb_name,"w+")
    outline = ''
    for x in hmhits:#extract seq from db for align
        outline += x+"\n"
        outline += hot_dic.get(x)+"\n"
    print(outline,file = subdb)
    os.system("cat "+seqname+" "+subdb_name+" > seq_hot.fasta")
    print("File seq_hot.fasta written!"+"\n"+"Alignment starts!")
    return subdb_name


# In[ ]:





# In[ ]:


def Align_analyze(seqname):
    ali_filename = seqname+".ali"
    os.system("/usr/local/bin/mafft"+"  --localpair  --maxiterate 16 --clustalout "+subdb_name+" > "+ali_filename)
    ali_file = open(ali_filename)
    next(ali_file)
    seq_dic = {}
    mark = ""
    for line in ali_file:
        if line.startswith(" "):
            mark = mark + line[16:].replace("\n", "")
        if line.startswith("\n"):
            continue
        else:
            head = line[0:15]
            #print(head)
            seq = line[16:].replace("\n", "")
            #print(seq)
            try:
                seq_dic[head] += seq
            except KeyError:
                seq_dic[head] = seq
    return seq_dic, mark


# In[ ]:


seq_dic = Ali_analyze("/Users/jsun/MHET/mhet_30.ali")[0]
mark = Ali_analyze("/Users/jsun/MHET/mhet_30.ali")[1]
#print(Ali_analyze("mhet_30.ali")[1])

import numpy as np
import pandas as pd
import math
import matplotlib 
target = "sp|A0A0K8P6T7|P"
def ali_array(seq_dic,mark,target):
    l = len(mark)
    lst = []
    for k in seq_dic:
        ls = list(seq_dic[k])
        lst.append(ls)
    df = pd.DataFrame(lst)
    t = seq_dic[target]
    target = []
    mark_lst = []
    cdict = {}
    n = 0
    s = 0
    mapdic = {}
    for x in t:
        if x == "-":
            s = s + 1
        else:
            target += x
            mark_lst += mark[n]
            c = mark[n]
            mapdic[n] = s
            n = n + 1
            
            if c == ".":
                cdict[n] = x+"_."
            if c == ":":
                cdict[n] = x+"_:"
            if c == "*":
                cdict[n] = x+"_*"
    return df,target,mark_lst,cdict,mapdic
df = ali_array(seq_dic,mark,target)[0]
t = ali_array(seq_dic,mark,target)[1]
m = ali_array(seq_dic,mark,target)[2]
#print(ali_array(seq_dic,mark,"sp|A0A0K8P6T7|P")[3])
#print(t,"\n",m)

cuttoff = 0.9
n = 0
print("position","WT","MUT")
for i in range(len(df.columns)):
    dic = df[i][0:19].value_counts().to_dict()
    x = df[i][20]
    if x == "-":
        continue
    else:
        n = n + 1
        for k in dic:
            p = int(dic[k])/19
            if p > cuttoff:
                if k == "-":
                    continue
                if k == x:
                    continue
                else:
                    print(n,x,k)

