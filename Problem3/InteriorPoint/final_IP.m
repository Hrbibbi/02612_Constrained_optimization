close all; clear all;

function x_sol = primal_dual_ipm(g, A, b, l, u, x, epsilon)
   
    % equality lagrange multiplier
    mu = 0; 

    %dual variables for the inequality constraints
    %lambdaL for the lower bound (0 <= x)
    %lambdaU for the upper bound (x <= Pd_max)
    lambdaL = ones(length(x), 1);
    lambdaU = ones(length(x), 1);

    %slack variables for the inequality consrtaints
    sL = x - l;
    sU = u - x;

    for k = 1:100
        %% precomputations
        SL = diag(sL);
        SU = diag(sU);
        LL = diag(lambdaL);
        LU = diag(lambdaU);

        %constant centering parameter
        sigma = 0.1;

        %% duality measure
        s = (sL'*lambdaL + sU'*lambdaU)/(2*length(x));

        %% compute the residuals
        %staionarity
        rL = g - A'*mu - lambdaL + lambdaU;

        %feasibility
        rA = A*x - b;

        %complementarity
        rCL = SL*LL*ones(length(sL),1) - sigma * s * ones(length(x),1);
        rCU = SU*LU*ones(length(sU),1) - sigma * s * ones(length(x),1);


        %% check convergence
        fprintf('Iteration %3i, |r_L| = %6.4e, |r_A| = %6.4e \n', k, norm(rL), norm(rA))
        if norm(rL) < epsilon && norm(rA) < epsilon && abs(s) < epsilon
            fprintf("Interior point: iteration%3i problem converged\n", k)
            x_sol = x;
            return;
        end
      
        %% newton step
        [n, m] = size(A);
        
        KKT = [zeros(m), -A',      -eye(m),     eye(m);
               A,        zeros(n), zeros(n,m),  zeros(n,m);
               LL,       zeros(m,n), SL,        zeros(m);
               LU,       zeros(m,n), zeros(m), SU];

        rhs = -[rL; rA; rCL; rCU];

        sol = KKT \ rhs;
        
        dx = sol(1:length(rL));
        dmu = sol(length(rL)+1 : length(rL)+length(rA));
        dlambdaL = sol(length(rL)+length(rA)+1 : length(rL)+length(rA)+length(rCL));
        dlambdaU = sol(length(rL)+length(rA)+length(rCL)+1:end);

        %% select alpha
        alpha = [-x./dx;-lambdaL./dlambdaL;-lambdaU./dlambdaU];
        num_zeros=alpha([dx<0;dlambdaL<0;dlambdaU<0]);
        if size(num_zeros,1)==0 
            alpha = 0;
        else
            alpha = min(num_zeros);
        end

        alpha = alpha*0.99;

        %% update
        x = x + alpha*dx;
        sL = x -l;
        sU = u - x;
        mu = mu + alpha*dmu;
        lambdaL = lambdaL + alpha*dlambdaL;
        lambdaU = lambdaU + alpha*dlambdaU;

    end
    fprintf("Did not converge\n");
    x_sol = x;
end

%% load variables for the problem
load("LP_Test.mat");
Cn = length(C);
Un = length(U);

g = [C;-U];
lenx = length(g);
A = [-ones(Cn,1);ones(Un,1)]';
b = 0;
u = [Pg_max; Pd_max];
l = [zeros(size(u))];
x0 = [(Un/Cn)*ones(Cn,1); ones(Un,1)];


x_sol = primal_dual_ipm(g, A, b, l, u, x0, 1e-8);

objective = g'*x_sol;
