import os
path = os.listdir(os.getcwd())
mp = os.getcwd()
plst = []
for p in path:
    if os.path.isdir(p):
        if len(p) == 2:
            plst.append(p)

for x in plst:
    os.chdir(x)
    #print(x)
    filename = "Dif_E-dimer_Repair.fxout"
    dif = open(filename)
    pdb_Dif = []
    for line in dif:
        if line.startswith("E-dimer"):
            lst = line.split("\t")
            t = (lst[0],lst[1])
            #print(t)
            pdb_Dif.append(t)
    file = "individual_list__"+x
    lines = []
    with open(file, "r") as f:
        for l in f:
            lines.append(l.replace("\n",""))
        f.close()
    n = 0
    while n < 500:
        try:
            v = float(pdb_Dif[n][1])
            #print(v)
            if v < -0.5:
                #print(x,pdb_Dif[n],lines[n])
                mut3 = lines[n].split(",")[0][-1]+lines[n].split(",")[1][-1]+lines[n].split(",")[2][-1]
                ddg = ("%.4f" % float(pdb_Dif[n][1]))
                print("mut"+mut3+".pdb"+"\t"+ddg)
                #os.system("cp "+pdb_Dif[n][0]+" "+mp+"/mut"+mut3+".pdb")
            n = n + 1
        except IndexError:
            os.chdir(mp)
            n = n + 1
    os.chdir(mp)
