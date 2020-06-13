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
