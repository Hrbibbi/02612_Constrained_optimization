clear
clc
load("QP_Test.mat")

%-----------------------------------
%           Quadprog part
%-----------------------------------

% A=[C';-C'];
% b=[du;-dl];
% options = optimoptions('quadprog', ...
%     'Display','iter', ...
%     'OptimalityTolerance',1e-4, ...
%     'ConstraintTolerance',1e-4 ...
% );
% tic
% [U,fval,exitflag,output,lambda]=quadprog(H,g,A,b,[],[],l,u,u,options);
% toc
% PlotSolutionQP(U)

%-----------------------------------
%           interior point
%-----------------------------------
tic
U=Interiorpointuppperlower(u,H,g,C,l,u,dl,du);
toc
minimizer=1/2*U'*H*U+g'*U;
PlotSolutionQP(U)
%-----------------------------------
%           Active set
%-----------------------------------

%U=activesetupperlower(H,g,C,l,u,dl,du);
%PlotSolutionQP(U)