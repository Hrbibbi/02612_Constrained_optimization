function x0 = feaseiblePointSimplex(A, b)
    % Find an initial feasible solution using the two-phase simplex method

    [m, n] = size(A);

    % Step 1: Try to solve the system A * x = b to check feasibility
    %disp('Checking feasibility of the system A * x = b...');
    try
        % Try to solve Ax = b directly (no artificial variables)
        x_initial = A \ b;
        if all(x_initial >= 0)
            % If the solution is non-negative, it's a feasible solution
            %disp('A feasible solution found directly without the need for artificial variables.');
            x0 = x_initial;  % Return this solution as the initial feasible point
            return;
        else
            %disp('Direct solution contains negative values, proceeding with Phase 1...');
        end
    catch
        %disp('No solution exists for A * x = b, proceeding with Phase 1...');
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
        %disp('Feasible solution found in Phase 1');
        x0 = x_phase1_sol(1:n);
    else
        error('No feasible solution found after Phase 1');
    end
end