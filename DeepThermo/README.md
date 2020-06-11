Background:
For a given protein sequence, there is a structure. Several mutations, within 10%, will not affect the overall structure but can drastically change many properties. For a typical enzyme of 300 amino acid (AA) long, exciting 300*19 single point mutations and 300*299*19 double mutations. To achieve a hyper-thermostable enzyme often requires about 10~15 mutations. Structure-based computation methods like FoldX, RosettaDDG, and AUACUS can predict stability change caused by mutations with modest accuracy. Computation based methods made achievement in protein thermostability design. It is not perfect but good enough compared to directed evolution.

For every mutation, software returns a score and if the score is within the cutoff, the mutation will be experimentally characterized for dTm. Here, in our previous work, 87 mutations were characterized, the dTm and the calculated score did not show any clear correlation and there is not any evidence or theory that the score and dTm are strongly correlated quantitatively, the cutoff is only qualitative measurement of dTm.

Aims:
Build a regression model to predict dTm using a score.
Create a classifier using more structure information to qualitatively predict dTm better.
Build a model for multi-site mutation prediction.
