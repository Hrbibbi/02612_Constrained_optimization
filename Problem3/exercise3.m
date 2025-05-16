load("LP_Test.mat")

g = [C; -U];
A = [-ones(15,1); ones(30,1)];
b = 0;
l = [zeros(45,1)];
u = [Pg_max; Pd_max];

f = g;
Aeq = A';
beq = b;
lb = l;
ub = u;

tic;
options = optimoptions('linprog', 'Display', 'iter');
[x_linprog, fval_linprog, exitflag, output, lambda] = ...
    linprog(f, [], [], Aeq, beq, lb, ub, options);
t_solve = toc;

Pg = x_linprog(1:15);
Pd = x_linprog(16:end);

disp('Total generation (MW):'), disp(sum(Pg))
disp('Total demand (MW):'), disp(sum(Pd))
disp('Market clearing price (€/MW):'), disp(lambda.eqlin)

sorted_C = sort(C);
sorted_U = sort(U,'descend');

cum_Pg = cumsum(Pg_max);
cum_Pd = cumsum(Pd_max);

figure;
hold on;
stairs([0; cum_Pg], [0; sorted_C], 'b', 'LineWidth', 2)
stairs([0; cum_Pd], [0; sorted_U], 'r', 'LineWidth', 2)
xlabel('Energy quantity (MW)')
ylabel('Price (€/MW)')
legend('Supply', 'Demand')
title('Market Clearing Supply-Demand Curve')
grid on;

fprintf('Solver stats:\n');
fprintf('Iterations: %d\n', output.iterations);
fprintf('CPU Time: %.4f seconds\n', t_solve);

print(gcf, 'supply_demand_plot', '-dpng', '-r300');  % 300 dpi for high quality



