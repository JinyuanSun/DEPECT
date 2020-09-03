#!/usr/bin/python3

import os
import sys

def cartddg_filter(ddgout):
    for line in ddgout:
        lst = line.split(":")
        if lst[2] == " WT":
            refdg = float(lst[3].split()[0])
        else:
            mutdg = float(lst[3].split()[0])
            ddg = mutdg - refdg
            mutname = lst[2].replace(" ","").split("_")[1]
            outlist.append(mutname+"\t"+str(round(ddg,2))+"\n")
    return outlist

def run(outlist):
    for mut_ddg in outlist:
        print(mut_ddg,end="")

if __name__ == '__main__':
    ddgout = open(sys.argv[1])
    outlist = []
    run(cartddg_filter(ddgout))
