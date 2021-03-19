# smiAAA
Stable Multi Input AAA
This repo contains code for generating the figures and data presented in "MULTI INPUT MINIMAX ADAPTIVE ANTOULAS-ANDERSON ALGORITHM
FOR RATIONAL APPROXIMATION WITH STABLE POLES".
*********
*Scripts*
*********
ch3_4_5_6.m : Runs the miAAA, smiAAA, miAAA-L, smiAAA-L and partial fraction conversion from Chapters 3,4,5,6.
Scripts Required:miaaa.m, smiaaa.m, przd.m, properrational.m, 
data:fdne.txt

holomorph.m : Runs the holomorphic function on the unit circle exmample from Chapter 5.
Scripts Requires:miaaa.m, circlemiaaa.m, properrational.m, przd.m

vfug_example.m : Runs the vector fitting users guide example from Chapter 8. Computers both the data in the table and the figure.
Scripts Required: symmetricsmiaa.m, properrational.m, przd.m
Toolboxes Required: RKFIT toolbox, Matrix Fitting Toolbox,
data:fdne.txt

iss_example.m : Runs the international space station example from Chapter 8. Computers both the data in the table and the figure.
Scripts Required: symmetricsmiaa.m, properrational.m, przd.m
Toolboxes Required: RKFIT toolbox, Matrix Fitting Toolbox,
data:iss.mat

mtf_example.m : Runs the Matrix Fitting Toolbox from Chapter 8. Computers both the data in the table and the figure.
Scripts Required: symmetricsmiaa.m, properrational.m, przd.m
Toolboxes Required: RKFIT toolbox, Matrix Fitting Toolbox,
data:ex2_Y.mat