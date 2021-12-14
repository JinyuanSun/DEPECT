#!/usr/bin/python3

ads = open("AutoDock41_scores.txt")
adslist = []
for line in ads:
    adslist.append(line.split()[-6:])

depectcafile = open("depect_enzde_geo.sc")
with open("depect_enzde_combine.sc", "w+") as newoutfile:
    i = 0
    for line in depectcafile:
        newoutfile.write(line.replace("\n", "," + ",".join(adslist[i])) + "\n")
        i = i + 1
