function [x_sol, objective, times] = InteriorPoint(g, A, b, l, u, x0, epsilon)
    warning('off', 'MATLAB:nearlySingularMatrix');
    %% convert to standard form
    %get the size of x
    n = length(x0);

    % make the slack variable vector for all inequalities x <= u
    s = ones(n, 1);

    %now we update the variables to be in standard form
    A = [A', zeros(1, n)
             eye(n), eye(n)];

    b = [b;u];
    g = [g; zeros(n, 1)];
    x = [x0;s];

    %the program is now in standard LP form
    % min g'x; A'x = b; x>=0;

    %% interior point algorithm

    lambda = ones(size(x));
    mu = zeros(size(b));

    %centering parameter
    sigma = 0.1;
    
    objective = g'*x;
    times = 0;

    for k = 1:100
        tic;
        %compute the duality measure
        s = (x'*lambda)/length(x);

        %compute residuals
        rL = g - A'*mu - lambda;
        rA = A*x - b;

        %check convergence
        %fprintf('Iteration %3i, |r_L| = %6.4e, |r_A| = %6.4e \n', k, norm(rL), norm(rA))
        if norm(rL) < epsilon && norm(rA) < epsilon && abs(s) < epsilon
            fprintf("Interior point: iteration%3i problem converged\n", k)
            break
        end

        %solve the newton step
        lenx = length(x);
        lenA = size(A, 1);

        diaglambda = diag(lambda);
        diagX = diag(x);

        KKT = [zeros(lenx),   -A',             -eye(lenx);
               A,          zeros(lenA),      zeros(lenA, lenx);
               diaglambda, zeros(lenx,lenA),       diagX];

        e = ones(lenx, 1);

        rhs = [-rL; -rA; -diagX*diaglambda*e+ sigma*s*e];

        sol = KKT \ rhs;

        dx = sol(1:length(rL));
        dmu = sol(length(rL)+1:length(rL)+length(rA));
        dlambda = sol(length(rL)+length(rA)+1:end);

        %compute alpha such that 
        %x + dx*alpha > 0 and lambda + dlambda*alpha > 0
        alpha = [-x./dx;-lambda./dlambda];
        num_zeros=alpha([dx<0;dlambda<0]);
        if size(num_zeros,1)==0 
            alpha = 0;
        else
            alpha = min(num_zeros);
        end

        alpha = alpha*0.99;

        x = x + dx*alpha;
        mu = mu + dmu*alpha;
        lambda = lambda + dlambda*alpha;

        times = [times, times(end)+toc];
        objective = [objective, g'*x];
    end
    times = times*1000;
    objective = objective*-1;
    x_sol = x(1:n);
    warning('on', 'MATLAB:nearlySingularMatrix');
end

