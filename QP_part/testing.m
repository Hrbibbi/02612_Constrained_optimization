
clear
clc

load('QP_Test.mat')
x=Interiorpointuppperlower(u,H,g,C,l,u,dl,du);
PlotSolutionQP(x)