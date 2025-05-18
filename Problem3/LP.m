close all; clear all;
%% add the function files to path
addpath("ActiveSet\")
addpath("InteriorPoint\")

%% load variables for the problem
load("LP_Test.mat");
Cn = length(C);
Un = length(U);

g = [C;-U];
A = [-ones(Cn,1);ones(Un,1)];
b = 0;
u = [Pg_max; Pd_max];
l = [zeros(size(u))];

%% convert to standard LP form
[g_std, A_std, b_std] = LPstandardForm(g, A, b, u);

%% exercise 3 (linprog)
tic;
options = optimoptions('linprog');
[x_linprog, fval_linprog, exitflag, output, lambda] = ...
    linprog(g, [], [], A', b, l, u, options);
time_linprog = toc;
iterations_linprog = output.iterations;

%% Interior point
x0 = [(Un/Cn)*ones(Cn,1);ones(Un,1)];
epsilon = 1e-8;
[x_IP, objective_IP, times_IP] = InteriorPoint(g,A,b,l,u, x0, epsilon);

%% Active set

x0 = feasiblePointLinprog(A_std, b_std);
[x_AS, objective_AS, times_AS] = Simplex(g_std, A_std, x0);

%% plotting

figure;
plot(times_IP, objective_IP, 'b-', 'LineWidth', 2); hold on;
plot(times_AS, objective_AS, 'r-', 'LineWidth', 2);

xlabel('Time (ms)');
ylabel('Social Welfare');
legend('Interior Point', 'Simplex', 'Location', 'southeast');
grid on;
set(gca, 'FontSize', 12);
saveas(gcf, 'welfare_plot.png');