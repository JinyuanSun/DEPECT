#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
import pandas as pd
import os
def _3D_coor21D(pdbfile):
    lst = []
    for line in pdbfile:
        if line.startswith("ATOM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            lst += [x,y,z]
    return lst


path = "/Users/jsun/RESEARCH/E-dimer/mutf" #文件夹目录
files = os.listdir(path)
for file in files:
    if file.endswith(".pdb"):
        pdbfile = open(file)
        lst = _3D_coor21D(pdbfile)
        #print(lst)
        try:
            df.loc[len(df)] = lst
        except NameError:
            df = pd.DataFrame([lst])
            
def eight(m):
    sp = "        "
    a = ("%.3f" % m)
    l = len(a)
    if l < 8:
        x = sp[0:8-l]+a
    else:
        x = a
    return x
f = open("E-dimer_Repair_1_67.pdb")
i = 0
newline = ''
of = open("ae.pdb","w+")
for line in f:
    if line.startswith("ATOM"):
        newline += line[0:30]+eight(df[i].mean())+eight(df[i+1].mean())+eight(df[i+2].mean())+"\n"
        i = i + 3
    else:
        newline += line
print(newline,file=of) 


# In[ ]:




