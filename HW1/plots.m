clc;clear;close all;

D = readtable("output.csv");

plot(D.Var1, D.Var2,'linewidth',2);