close all; clear all;
%% add the function files to path
addpath("ActiveSet\")
addpath("InteriorPoint\")

%% load variables for the problem
load("LP_Test.mat");
Cn = length(C);
Un = length(U);

g = [C;-U];
lenx = length(g);
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
time_linprog = toc*1000;
iterations_linprog = output.iterations;

%% Interior point
x0 = [(Un/Cn)*ones(Cn,1);ones(Un,1)];
epsilon = 1e-8;
[x_IP, objective_IP, times_IP, duals_IP, primals_IP, complementarities_IP] = InteriorPoint(g_std,A_std,b_std, x0, epsilon);

%% Active set
tic;
x0 = feaseiblePointSimplex(A_std, b_std);
feasible_time = toc*1000;
[x_AS, objective_AS, times_AS, duals_AS, primals_AS, complementarities_AS] = Simplex(g_std, A_std, b_std, x0);
x_AS = x_AS(1:lenx);
times_AS = times_AS + ones(size(times_AS))*feasible_time;

%% Time plot
figure;
plot(times_IP, objective_IP, 'b-', 'LineWidth', 2); hold on;
plot(times_AS, objective_AS, 'r-', 'LineWidth', 2);

xlabel('Time (ms)');
ylabel('Social Welfare');
legend('Interior Point', 'Simplex', 'Location', 'southeast');
grid on;
set(gca, 'FontSize', 12);
saveas(gcf, 'welfare_plot.png');

%% IP residuals plot

figure;
iters = 1:length(primals_IP);

%primal
subplot(2,3,1);
plot(iters, primals_IP, 'r-', 'LineWidth', 2);
title('Interior Point - Primal residual');
xlabel('Iteration'); ylabel('r_A'); grid on;

%dual
subplot(2,3,2);
plot(iters, duals_IP, 'b-', 'LineWidth', 2);
title('Interior Point - Dual residual');
xlabel('Iteration'); ylabel('r_L'); grid on;

%complimentarity
subplot(2,3,3);
plot(iters, complementarities_IP, 'k-', 'LineWidth', 2);
title('Interior Point - Complementarity');
xlabel('Iteration'); ylabel('r_{SZ}'); grid on;

%% Simplex residuals plot

%primal
subplot(2,3,4);
plot(iters, primals_AS, 'r-', 'LineWidth', 2);
title('Simplex - Primal residual');
xlabel('Iteration'); ylabel('r_A'); grid on;

%dual
subplot(2,3,5);
plot(iters, duals_AS, 'b-', 'LineWidth', 2);
title('Simplex - Dual residual');
xlabel('Iteration'); ylabel('r_L'); grid on;

%complimentarity
subplot(2,3,6);
plot(iters, complementarities_AS, 'k-', 'LineWidth', 2);
title('Simplex - Complementarity');
xlabel('Iteration'); ylabel('r_{SZ}'); grid on;

%% IP & Simplex - Objective plot

iters_IP = 1:length(objective_IP);
iters_AS = 1:length(objective_AS);

figure;

%interior point
subplot(1,2,1);
ytickformat('%.1e')
plot(iters_IP, objective_IP, 'm-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Social welfare (objective value)');
title('Interior point');
grid on;
set(gca, 'FontSize', 12);

%simplex
subplot(1,2,2);
plot(iters_AS, objective_AS, 'gint-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Social welfare (objective)');
title('Simplex');
grid on;
set(gca, 'FontSize', 12);