function x0 = feasiblePointLinprog(A, b)
    % This function finds an initial feasible point for the simplex algorithm using linprog

    % A: Constraint matrix (m x n)
    % b: Right-hand side vector (m x 1)

    [m, n] = size(A); 

    % Objective function for the auxiliary linear program
    f = [zeros(n,1); 1];  

    % Construct the augmented matrix A_lp for the linear program
    A_lp = [ A, -ones(m,1);   % First block: A with slack variable of +1
             -A, -ones(m,1) ]; % Second block: -A with slack variable of -1

    % Construct the right-hand side for the LP
    b_lp = [ b; -b ];

    % Set the lower bounds for the LP variables
    lb = [zeros(n,1); 0];

    % Solve the linear program using MATLAB's linprog function
    z = linprog(f, A_lp, b_lp, [], [], lb, []);  % Solve LP

    % Extract the original feasible solution from the solution vector z
    x0 = z(1:n);  % Extract the original variable values from z

    % Extract the optimal value of the slack variable
    t_opt = z(end);  % The last value in z corresponds to the slack variable's value

    % If the optimal value of the slack variable is positive, the system is not feasible
    if t_opt > 0
        error("not feasible")  % If t_opt > 0, there is no feasible solution
    end
end
