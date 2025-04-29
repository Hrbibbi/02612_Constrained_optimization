function [x_opt, W_opt] = active_set_QP(G, c, A, b, x0, W0, max_iter)
    % Inputs:
    % G: n x n positive semi-definite matrix
    % c: n x 1 vector for the linear terms of the objective
    % A: m x n matrix for the inequality constraint matrix
    % b: m x 1 vector for the inequality constraint bounds
    % x0: Initial feasible point (n x 1)
    % W0: Initial active set of inequality indices (set of active constraints)
    % max_iter: Maximum number of iterations
    
    % Initialize variables
    x_k = x0;
    W_k = W0;
    
    for k = 1:max_iter
        % Step 1: Solve the QP subproblem for the step pk
        [p_k, lambda_k] = solve_qp_subproblem(G, c, A, x_k, W_k);
        
        if norm(p_k) < 1e-6  % If pk is very close to zero
            % Step 2: Compute Lagrange multipliers and check the KKT conditions
            [lambda, W_opt] = compute_lagrange_multipliers(G, c, A, x_k, W_k);
            
            % Check if we have a solution satisfying the KKT conditions
            if all(lambda >= 0)
                x_opt = x_k;
                return;
            else
                % Step 3: Remove the most negative multiplier and update the working set
                [min_lambda, j] = min(lambda(W_k));
                W_k = setdiff(W_k, j);
            end
        else
            % Step 4: Compute the step length
            alpha_k = compute_step_length(A, b, x_k, p_k, W_k);
            
            % Step 5: Update the solution
            x_k = x_k + alpha_k * p_k;
            
            % Step 6: Check for blocking constraints and update the working set
            if alpha_k < 1
                new_constraint = find_new_constraint(A, b, x_k, p_k, W_k);
                W_k = union(W_k, new_constraint);
            end
        end
    end
    
    x_opt = x_k;
    W_opt = W_k;
end

function [p_k, lambda_k] = solve_qp_subproblem(G, c, A, x_k, W_k)
    % Solves the equality-constrained QP subproblem at iteration k
    % min_p  0.5 * p' * G * p + g_k' * p
    % subject to A(W_k)' * p = 0 (W_k constraints are active)
    
    % Gradient of the objective function at x_k
    g_k = G * x_k + c;
    
    % Submatrix of A corresponding to the active set W_k
    A_Wk = A(W_k, :);
    
    % Solve the equality-constrained QP
    p_k = (G(W_k, W_k)) \ (-g_k(W_k));
    lambda_k = A_Wk' * p_k;  % Calculate multipliers
end

function [lambda, W_opt] = compute_lagrange_multipliers(G, c, A, x_k, W_k)
    % Computes the Lagrange multipliers for the current working set
    g_k = G * x_k + c;
    A_Wk = A(W_k, :);
    
    % Solve for multipliers that satisfy the KKT conditions
    lambda = (A_Wk' * A_Wk) \ A_Wk' * g_k;
    W_opt = W_k;
end

function alpha_k = compute_step_length(A, b, x_k, p_k, W_k)
    % Computes the step length alpha_k based on the blocking constraints
    alpha_k = 1;
    for i = 1:length(W_k)
        if A(i, :) * p_k < 0
            alpha_k = min(alpha_k, (b(i) - A(i, :) * x_k) / (A(i, :) * p_k));
        end
    end
end

function new_constraint = find_new_constraint(A, b, x_k, p_k, W_k)
    % Find a new blocking constraint to add to the working set
    for i = 1:size(A, 1)
        if A(i, :) * p_k < 0 && ~ismember(i, W_k)
            new_constraint = i;
            return;
        end
    end
end
