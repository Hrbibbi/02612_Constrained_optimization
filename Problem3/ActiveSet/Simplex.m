function [x_sol, objective, times] = Simplex(g, A, x)
    
    %create set N of all i where xi = 0
    %and set B with all the other xi (xi > 0)
    Fullset = 1:length(x);
    Nset = find(x==0);
    Bset = setdiff(Fullset,Nset);

    objective = g'*x;
    times = 0;

    for k = 0:100
        tic;

        xB = x(Bset);
        xN = x(Nset);
        B = A(:,Bset);
        N = A(:,Nset);
        mu = B'\g(Bset);
        lambdaN = g(Nset) - N'*mu;

        if all(lambdaN>=0)
            fprintf("Simplex: iteration %i problem converged\n", k)
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

        times = [times, times(end)+toc];
        objective = [objective, g'*x];
    end
    times = times*1000;
    objective = objective*-1;
    x_sol = x;
end