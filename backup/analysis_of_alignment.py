#!/usr/bin/env python3
# coding: utf-8

def Ali_analyze(ali_filename):
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