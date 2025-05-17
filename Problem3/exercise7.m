close all; clear all;

function [g_std, A_std, b_std] = LPstandardForm(g, A, b, u)
    %convert the lp problem into standard form

    [m, n] = size(A);

    %now we update the variables to be in standard form
    g_std = [g; zeros(m, 1)];

    A_std = [A', zeros(n, m)
             eye(m), eye(m)];

    b_std = [b;u];
end

% function x0 = feasiblePoint(A, b)
%     [m, n] = size(A);
%     f = [zeros(n,1); 1];                        % objective: min t
%     A_lp = [ A, -ones(m,1);                   % Ax - t <= b
%             -A, -ones(m,1) ];                 % -Ax - t <= -b
%     b_lp = [ b; -b ];
% 
%     lb = [zeros(n,1); 0];                      % x ≥ 0, t ≥ 0
% 
%     z = linprog(f, A_lp, b_lp, [], [], lb, []);
%     x0 = z(1:n);
%     t_opt = z(end);
% 
%     if t_opt > 0
%         error("not feasible")
%     end
% end

function x0 = feasiblePoint(A, b)

    [m, n] = size(A);


    %make x0 = [x; t; s1; s2]
    x = zeros(n, 1);

    t = max(abs(b));

    s1 = t*ones(m, 1) - b;
    s2 = t*ones(m, 1) + b;

    x0 = [x; t; s1; s2];

    %create g so that it minimizes t only
    g = [zeros(n, 1); 1; zeros(2*m, 1)];

    %create A to satisfy the new equality conditions
    % Ax + t >= +-b     Ax + t - s = +-b
    A = [A, ones(m, 1), -eye(m), zeros(m);
        -A, ones(m, 1), zeros(m), -eye(m)];

    %find the feasible point using simplex
    x0 = simplex(g, A, x0);
    t = x0(n+1);
    x0 = x0(1:n);

    if t > 1e-6
        error("Infeasible")
    end
end

function x_sol = simplex(g, A, x)
    
    %create set N of all i where xi = 0
    %and set B with all the other xi (xi > 0)
    Fullset = 1:length(x);
    Nset = find(x==0);
    Bset = setdiff(Fullset,Nset);

    for k = 0:100
        xB = x(Bset);
        xN = x(Nset);
        B = A(:,Bset);
        N = A(:,Nset);
        mu = B'\g(Bset);
        lambdaN = g(Nset) - N'*mu;

        fprintf("t = %.4f\n", x(90+1));

        if all(lambdaN>=0)
            fprintf("iteration %i: problem converged\n", k)
            break
        end
        
        s = find(lambdaN<0, 1);
        h = B\N(:,s);
        idx = 1:length(xB);
        idx = idx(h>0);
        sol = xB(h>0)./h(h>0);
        [~, J] = min(sol);
        J = idx(J);
        if J == []
            error("Unbound problem, no solution\n")
            return;
        end

        alpha = xB(J)./h(J);
        xB = xB - alpha*h;
        xB(J) = 0;
        xN(s) = alpha;
        x(Bset) = xB;
        x(Nset) = xN;
        temp = Bset(J);
        Bset(J) = Nset(s);
        Nset(s) = temp;    
        Bset = sort(Bset);
        Nset = sort(Nset);
    end
    x_sol = x;
end

%load variables for the problem
load("LP_Test.mat");

Cn = length(C);
Un = length(U);

g = [C;-U];
A = [-ones(Cn,1);ones(Un,1)];
b = 0;
u = [Pg_max; Pd_max];
l = [zeros(size(u))];


%transform problem into standard form
[g, A, b] = LPstandardForm(g, A, b, u);

%get feasible point (using simplex)
fprintf("phase 1\n")
x0 = feasiblePoint(A, b);


%get an optimal point (using simplex)
fprintf("phase 2\n")
x = simplex(g, A, x0);