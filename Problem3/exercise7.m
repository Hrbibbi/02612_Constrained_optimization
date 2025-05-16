close all; clear all;


function x_sol = simplex(g, A, b, l, u, x0)
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

    %% simplex algorithm
    %create set N of all i where xi = 0
    %and set B with all the other xi (xi > 0)
    Fullset = 1:length(x);
    Nset = find(x==0);
    Bset = setdiff(Fullset,Nset);

    %remake the A matrix as well
    Bmat = A(:, Bset);
    Nmat = A(:, Nset);

    A = [Bmat, Nmat];

    %calcualte xb
    xb = Bmat \ b;
    xn = zeros(length(Nset), 1);

    %reorder x to be [xb;xn]
    x = zeros(length(x), 1);
    x(Bset) = xb;
    x(Nset) = xn;

    %set lambda accordingly
    %for all xb lambda = 0 and for all xn lambda >= 0
    lambda = zeros(length(x), 1);
    lambda(Nset) = 1;
    lambda(Bset) = 0;

    for k = 0:100
        xB = x(Bset);
        xN = x(Nset);
        B = A(:,Bset);
        N = A(:,Nset);
        mu = B'\g(Bset);
        lambdaN = g(Nset) - N'*mu;

        if all(lambdaN>=0)
            fprintf("Problem converged\n")
            break
        end
        
        s = find(lambdaN<0, 1);
        h = B\N(:,s);
        idx = 1:length(xB);
        idx = idx(h>0);
        sol = xB(h>0)./h(h>0);
        min_sol = min(sol);
        J = min_sol(2);
        J = index(J);
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
x0 = [(Un/Cn)*ones(Cn,1);ones(Un,1)];

x = simplex(g, A, b, l, u, x0);