Background:
For a given protein sequence, there is a structure. Several mutations, within 10%, will not affect the overall structure but can drastically change many properties. For a typical enzyme of 300 amino acid (AA) long, exciting 300X19 single point mutations and 300X299X19 double mutations. To achieve a hyper-thermostable enzyme often requires about 10~15 mutations. Structure-based computation methods like FoldX, RosettaDDG, and AUACUS can predict stability change caused by mutations with modest accuracy. Computation based methods made achievement in protein thermostability design. It is not perfect but good enough compared to directed evolution.

For every mutation, software returns a score and if the score is within the cutoff, the mutation will be experimentally characterized for dTm. Here, in our previous work, 87 mutations were characterized, the dTm and the calculated score did not show any clear correlation and there is not any evidence or theory that the score and dTm are strongly correlated quantitatively, the cutoff is only qualitative measurement of dTm.

Aims:
Build a regression model to predict dTm using a score.
Create a classifier using more structure information to qualitatively predict dTm better.
Build a model for multi-site mutation prediction.


-------  
for h5li
Design:  
旋转蛋白质的结构(把突变氨基酸附近（6Å）范围内的氨基酸拿出来，或者用整个结构)来实现数据增强；希望达到弱化背景的某些特点，因为结构的信息相当丰富。 
对于同一个突变不同软件可能会预测出不同的局部结构，也会给出不一样的打分，利用结构的三维坐标、氨基酸和原子类型、打分和实验值来训练神经网络，提高预测的准确率。
hydropathy_index = {"R":-2.5,"K":-1.5,"D":-0.9,"Q":-0.85,"N":-0.78,"E":-0.74,"H":-0.4,"S":-0.18,"T":-0.05,
                        "P":0.12,"Y":0.26,"C":0.29,"G":0.48,"A":0.62,"M":0.64,"W":0.81,"L":1.1,"V":1.1,
                        "F":1.2,"I":1.4}
occurance = {"A":8.76,"R":5.78,"N":3.93,"D":5.49,"C":1.38,"Q":3.9,"E":6.32,"G":7.03,"H":2.26,"I":5.49,"L":9.68,
                 "K":5.19,"M":2.32,"F":3.87,"P":5.02,"S":7.14,"T":5.53,"W":1.25,"Y":2.91,"V":6.73}
charge = {"A":0,"D":-1,"E":-1,"H":+0.1,"C":0,"Y":0,"K":+1,"R":+1,"N":0,"Q":0,"G":0,"I":0,"L":0,"M":0,"F":0,"P":0,
              "S":0,"T":0,"W":0,"V":0}
