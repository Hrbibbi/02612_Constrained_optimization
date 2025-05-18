function x0 = feasiblePointSimplex(A, b)
    % Find an initial feasible solution using the two-phase simplex method

    [m, n] = size(A);

    % Step 1: Try to solve the system A * x = b to check feasibility
    disp('Checking feasibility of the system A * x = b...');
    try
        % Try to solve Ax = b directly (no artificial variables)
        x_initial = A \ b;
        if all(x_initial >= 0)
            % If the solution is non-negative, it's a feasible solution
            disp('A feasible solution found directly without the need for artificial variables.');
            x0 = x_initial;  % Return this solution as the initial feasible point
            return;
        else
            disp('Direct solution contains negative values, proceeding with Phase 1...');
        end
    catch
        disp('No solution exists for A * x = b, proceeding with Phase 1...');
    end

    % Step 2: Phase 1: Construct artificial variables
    A_phase1 = [A, eye(m)];
    g_phase1 = [zeros(n, 1); ones(m, 1)];  % Minimize the sum of artificial variables
    b_phase1 = b;

    % Initial basic solution will be the artificial variables
    x_phase1 = [zeros(n, 1); b_phase1];  % Initially, all x are zero, but artificial variables will be non-zero
    
    % Use Simplex to solve for Phase 1
    [x_phase1_sol, ~, ~] = Simplex(g_phase1, A_phase1, x_phase1);
    
    % Check if the objective function is 0, meaning we found a feasible point
    if all(x_phase1_sol(n+1:end) == 0)  % If no artificial variables are in the basis, we have a feasible solution
        % Extract the basic feasible solution (x_phase1_sol corresponds to the original variables)
        disp('Feasible solution found in Phase 1');
        x0 = x_phase1_sol(1:n);
    else
        error('No feasible solution found after Phase 1');
    end
end


function [g_std, A_std, b_std] = LPstandardForm(g, A, b, u)
    % Convert the LP problem into standard form

    [m, n] = size(A);

    % Now we update the variables to be in standard form
    g_std = [g; zeros(m, 1)];

    A_std = [A', zeros(n, m);
             eye(m), eye(m)];

    b_std = [b; u];
end

function [x_sol, objective, times] = Simplex(g, A, x)
    % The primal Simplex algorithm for solving linear programs

    % Create set N of all i where xi = 0
    % and set B with all the other xi (xi > 0)
    Fullset = 1:length(x);
    Nset = find(x == 0);
    Bset = setdiff(Fullset, Nset);

    objective = g' * x;
    times = 0;

    for k = 0:100
        tic;

        xB = x(Bset);
        xN = x(Nset);
        B = A(:, Bset);
        N = A(:, Nset);
        mu = B' \ g(Bset);
        lambdaN = g(Nset) - N' * mu;

        if all(lambdaN >= 0)
            fprintf("Simplex: iteration %i problem converged\n", k);
            break;
        end
        
        s = find(lambdaN < 0, 1);
        h = B \ N(:, s);
        idx = 1:length(xB);
        idx = idx(h > 0);
        sol = xB(h > 0) ./ h(h > 0);
        [~, J] = min(sol);
        J = idx(J);
        if isempty(J)
            error("Unbound problem, no solution\n");
            return;
        end

        alpha = xB(J) / h(J);
        xB = xB - alpha * h;
        xB(J) = 0;
        xN(s) = alpha;
        x(Bset) = xB;
        x(Nset) = xN;
        temp = Bset(J);
        Bset(J) = Nset(s);
        Nset(s) = temp;    
        Bset = sort(Bset);
        Nset = sort(Nset);

        times = [times, times(end) + toc];
        objective = [objective, g' * x];
    end
    times = times * 1000;
    objective = objective * -1;
    x_sol = x;
end

% Main script to load the problem data and solve using Simplex
load("LP_Test.mat");
Cn = length(C);
Un = length(U);

g = [C; -U];
A = [-ones(Cn, 1); ones(Un, 1)];
b = 0;
u = [Pg_max; Pd_max];
l = [zeros(size(u))];

%% Convert to standard LP form
[g_std, A_std, b_std] = LPstandardForm(g, A, b, u);

% Find an initial feasible point using the Two-Phase Simplex method
x0 = feasiblePointSimplex(A_std, b_std);

% Solve the LP using the Primal Simplex algorithm
[x_AS, objective_AS, times_AS] = Simplex(g_std, A_std, x0);

% Display the solution and other details
disp('Optimal solution:');
disp(x_AS);
disp('Objective value:');
disp(objective_AS(end));
disp('Computation time (ms):');
disp(times_AS);
