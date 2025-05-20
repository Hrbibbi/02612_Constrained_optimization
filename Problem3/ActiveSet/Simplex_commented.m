function [x_sol, objective, times] = Simplex(g, A, x)
    % g: cost vector
    % A: equality constraint matrix
    % x: initial guess for the solution

    % Create the set of basic and non-basic variables
    Fullset = 1:length(x);
    Nset = find(x==0);
    Bset = setdiff(Fullset,Nset);

    objective = g'*x;
    times = 0;

    for k = 0:100  
        tic;  

        % Extract variables and matrices corresponding to B and N sets
        xB = x(Bset);   
        xN = x(Nset);   
        B = A(:,Bset);  
        N = A(:,Nset);  
      
        % Compute the dual variables from the B matrix and g vector
        mu = B'\g(Bset);
        
        % Compute reduced costs for non-basic variables
        lambdaN = g(Nset) - N'*mu;

        % Check if the problem has converged
        if all(lambdaN >= 0)
            fprintf("Simplex: iteration %i problem converged\n", k)
            break
        end
        
        % Select entering variable with the most negative reduced cost
        s = find(lambdaN < 0, 1);
        
        % Compute direction vector for the pivot operation
        h = B\N(:,s);
        
        % Find the leaving variable based on the minimum ratio test
        idx = 1:length(xB);
        idx = idx(h > 0);  
        sol = xB(h > 0)./h(h > 0);  
        [~, J] = min(sol);
        J = idx(J);

        % If no leaving variable is found, the problem is unbounded
        if J == []
            error("Unbound problem, no solution\n")
            return;
        end

        % Compute the step size for the pivot operation
        alpha = xB(J) / h(J);

        % Update the basic and non-basic variables
        xB = xB - alpha*h;  
        xB(J) = 0;
        xN(s) = alpha;

        % Update the solution vector with the new B and N values
        x(Bset) = xB;
        x(Nset) = xN;

        % Swap the entering and leaving variables between Bset and Nset
        temp = Bset(J);
        Bset(J) = Nset(s);
        Nset(s) = temp;    

        Bset = sort(Bset);
        Nset = sort(Nset);

        % Record the time taken for this iteration and the objective value
        times = [times, times(end) + toc];
        objective = [objective, g'*x];
    end

    times = times * 1000;
    
    % Reverse the sign of the objective (since we are minimizing)
    objective = objective * -1;
    
    % Return the final solution, objective values, and computation times
    x_sol = x;
end
