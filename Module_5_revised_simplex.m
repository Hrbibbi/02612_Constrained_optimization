clear all; close all;

function [x] = revisedSimplex(A,b,g,x)
    Fullset = 1:length(x);
    Nset = find(x==0);
    Bset = setdiff(Fullset,Nset);
    itMax = 100;
    for it = 1:itMax
        display([x(1)-x(2),x(3)-x(4)])
        xB = x(Bset);
        xN = x(Nset);
        B = A(:,Bset);
        N = A(:,Nset);
        mu = B'\g(Bset);
        lambdaN = g(Nset)-N'*mu;
        if all(lambdaN>=0)
            fprintf('Problem converged')
            break
        else
            s = (find(lambdaN<0,1));
            h = B\N(:,s);
            index = 1:length(xB);index = index(h>0);
            temp = xB(h>0)./h(h>0);
            [~,J] = min(temp);
            J = index(J);
            if J == [];
                error("Unbound problem, No solution")
            else
                alpha = xB(J)./h(J);
                xB = xB-alpha*h;
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
        end
    end
end

A = [1 -1 0  0 -1 0  0 0;
     1 -1 0  0  0 1  0 0;
     0 0  1 -1  0 0 -1 0;
     0 0  1 -1  0 0  0 1];

b = [-1 1 -1 1]';

g = [1 -1 -1 1 0 0 0 0]';

x0 = [1 0 0 1 2 0 0 2]';

[x] = revisedSimplex(A,b,g,x0);
