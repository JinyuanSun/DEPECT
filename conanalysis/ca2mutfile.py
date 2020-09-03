#!/usr/bin/python3


import os
import argparse

def ca2mutfilelist(depectcafile):
    cafile = open(depectcafile)
    i = 0
    mutfilelist = []
    for line in cafile:
        if line.startswith("#"):
            continue
        else:
            lst = line.split("\t")
            mut_unit = lst[0]+" "+lst[1]+" "+lst[2]
            i = i + 1
            mutfilelist.append(mut_unit)
    return mutfilelist


def mutfilelist2rosettamutfile(mutfilelist,mutfilename):
    #mutfilename = depectcafile.split(".")[0]+".mutfile"
    mutfile = open(mutfilename,"w+")
    with open(mutfilename,"a+") as mutfile:
        mutfile.write("total "+str(len(mutfilelist))+"\n")
        for mutunit in mutfilelist:
            mutfile.write("1\n"+mutunit+"\n")
        mutfile.close()


def chunk(lst,nt,prefix):
    outfilenamelist = []
    t = len(lst)//(int(nt)-1)
    n = 0
    for i in range(0,len(lst),t):
        b=lst[i:i+t]
        n = n + 1
        mutfilelist2rosettamutfile(b,prefix+"_"+str(n)+".mutfile")
        c = prefix+"_"+str(n)+".mutfile"
        outfilenamelist.append(c)
    return outfilenamelist




def get_parser():
    parser = argparse.ArgumentParser(description='convert depect consensus analysis output .ca file to rosetta mutfile')
    #parser.parse_args()
    parser.add_argument("-in", "--input", help="depect consensus analysis output .ca file")
    parser.add_argument("-nt", "--number_splited", help="number of output mutfile, default is 1",default=1)
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    depectcafile = args.input
    nt = int(args.number_splited)
    mutfilelist = ca2mutfilelist(depectcafile)
    #print(nt)
    outfilename = prefix = depectcafile.split(".")[0]
    if nt == 1:
        mutfilelist2rosettamutfile(mutfilelist,outfilename)
    if nt > 1:
        chunk(mutfilelist,nt,prefix)
