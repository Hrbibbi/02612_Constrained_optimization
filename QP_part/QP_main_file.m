clear
clc
close all
load("QP_Test.mat")

%-----------------------------------
%           Quadprog part
%-----------------------------------

A=[C';-C'];
b=[du;-dl];
options = optimoptions('quadprog', ...
    'Display','iter', ...
    'OptimalityTolerance',1e-4, ...
    'ConstraintTolerance',1e-4 ...
);
tic
[U_q,fval,exitflag,output,lambda]=quadprog(H,g,A,b,[],[],l,u,u,options);
toc
%PlotSolutionQP(U)
%plot(lambda)
lam_C=lambda.ineqlin(1:200)
lam_Cneg=lambda.ineqlin(201:end)
lam_l=lambda.lower
lam_u=lambda.upper

%-----------------------------------
%           interior point
%-----------------------------------

tic
[U_int,yl,yu,zl,zu]=Interiorpointuppperlower(u,H,g,C,l,u,dl,du);
toc

% PlotSolutionQP(U)
%-----------------------------------
%           Active set
%-----------------------------------

%PlotSolutionQP(U)

[x_act, lam_l_act, lam_u_act, lam_Cneg_act, lam_C_act] = activesetupperlower(H, g, C, l, u, dl, du);

figure;

% Row 1: quadprog
subplot(3,5,1);  plot(U_q, 'k');       title('quadprog: solution');          grid on;
subplot(3,5,2);  plot(lam_l, 'g');   title('\lambda_{lower}');            grid on;
subplot(3,5,3);  plot(lam_u, 'm');   title('\lambda_{upper}');            grid on;
subplot(3,5,4);  plot(lam_Cneg, 'b');title('\lambda_{C} (C^T x \geq dl)');grid on;
subplot(3,5,5);  plot(lam_C, 'r');   title('\lambda_{Cneg} (C^T x \leq du)'); grid on;

% Row 2: interior point
subplot(3,5,6);  plot(U_int, 'k');       title('int.point: solution');         grid on;
subplot(3,5,7);  plot(yl, 'g');      title('\lambda_{lower}');            grid on;
subplot(3,5,8);  plot(yu, 'm');      title('\lambda_{upper}');            grid on;
subplot(3,5,9);  plot(zl, 'b');      title('\lambda_{C} (C^T x \geq dl)');grid on;
subplot(3,5,10); plot(zu, 'r');      title('\lambda_{Cneg} (C^T x \leq du)'); grid on;

% Row 3: active set
subplot(3,5,11); plot(x_act, 'k');          title('activeset: solution');       grid on;
subplot(3,5,12); plot(lam_l_act, 'g');      title('\lambda_{lower}');           grid on;
subplot(3,5,13); plot(lam_u_act, 'm');      title('\lambda_{upper}');           grid on;
subplot(3,5,14); plot(lam_Cneg_act, 'b');   title('\lambda_{C} (C^T x \geq dl)');grid on;
subplot(3,5,15); plot(lam_C_act, 'r');      title('\lambda_{Cneg} (C^T x \leq du)'); grid on;

sgtitle('Comparison of Variables: quadprog vs Interior Point vs Active Set');