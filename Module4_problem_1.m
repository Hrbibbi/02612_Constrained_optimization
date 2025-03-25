clear all; close all;

function [H,g,A,b] = generateProblem(n, uBar, d0)
    H = speye(n+1);
    Ai = [1 1 n n n 2:n-1 2:n-1];
    Aj = [1 n n-1 n n+1 1:n-2 2:n-1];
    Ak = [-1; 1; 1; -1; -1;ones(n-2,1); repelem(-1,n-2,1)]; 
    A = sparse(Ai,Aj,Ak);
    g = repelem(uBar,n+1,1);
    b = zeros(n,1);
    b(1) = -d0;
end

function [S, F] = KKTequation(H,g,A,b)
    S = [H ,-A';-A,sparse(size(b,1),size(b,1))];
    F = -[g;b];
end

function [x, lambda] = solveNullSpace(H,g,A,b);
    [Q,Rbar] = qr(A');
    m1 = size(Rbar,2);
    Q1 = Q(:,1:m1);
    Q2 = Q(:,m1+1:size(A,1));
    R = Rbar(1:m1,1:m1);
    xy = R\b;
    xz = (Q2'*H*Q2)\(-Q2'*(H*Q1*xy+g));
    x = Q1*xy+Q2*xz;
    lambda = R\(Q1'*(H*x+g));
end

function [x, lambda] = solveRangeSpace(H,g,A,b);
    %% assumption that H is identity:
    Ha = A*A';
    lambda = Ha\(b+A*g);
    x = H\(A'*lambda-g);
end

uBar = 0.2; 
d0 = 1;
ns = floor(logspace(1,4,10));
for i =1:length(ns);
    [H,g,A,b] = generateProblem(ns(i), uBar, d0);
    Hfull = full(H); Afull = full(A);
    [S,f] = KKTequation(H,g,A,b);
    Sfull = full(S);

% LU factorization
tic
[L,U,P] = lu(S);
xLU = U\(L\(P*f));
tLU(i) = toc;

%LDL factorization
tic
[L,D,P] = ldl(S);
xLDL = P*(L'\(D\(L\(P'*f))));
tLDL(i) = toc;

%Range-Space
tic
[x, lambda] = solveRangeSpace(H,g,A,b);
tRangeSpace(i) = toc;

% LU factorization (full)
tic
[L,U,P] = lu(Sfull);
xLU = U\(L\(P*f));
tLUfull(i) = toc;

%LDL factorization (full)
tic
[L,D,P] = ldl(Sfull);
xLDL = P*(L'\(D\(L\(P'*f))));
tLDLfull(i) = toc;

%Range-Space  (full)
tic
[x, lambda] = solveRangeSpace(Hfull,g,Afull,b);
tRangeSpacefull(i) = toc;

end
%Null space
%[x, lambda] = solveNullSpace(H,g,A,b);

figure()
loglog(ns,tLU,':ok','DisplayName','LU-factorization (sparse)')
hold on;
loglog(ns,tLDL,':ob','DisplayName','LDL-factorization (sparse)')
loglog(ns,tRangeSpace,':or','DisplayName','Range-Space (sparse)')

loglog(ns,tLUfull,'--ok','DisplayName','LU-factorization (full)')
loglog(ns,tLDLfull,'--ob','DisplayName','LDL-factorization (full)')
loglog(ns,tRangeSpacefull,'--or','DisplayName','Range-Space (full)')
legend();
xlabel("Problem size n")
ylabel("Computation time [s]")