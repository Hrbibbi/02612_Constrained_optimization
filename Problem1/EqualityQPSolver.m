function [x,lambda] = EqualityQPSolver(H,g,A,b,solver)
% Equality Quadratic Problem Solver
    %   min     1/2*x'*H*x+g'*x
    %    x
    %   s.t.    A'x=b
    %
    % [x,lambda]=EqualityQPSolver(H,g,A,b,solver)
    %   Inputs:
    %       H           Matrix
    %       g           Vector
    %       A           Matrix
    %       b           Vector
    %       solver      solver type
    %
    %   Outputs: optimized x and langrange multipliers lambda:
    
    if strcmpi(solver,"LUdense")
        [x,lambda] = EqualityQPSolverLUdense(H,g,A,b);
    elseif strcmpi(solver,"LUsparse")
        [x,lambda] = EqualityQPSolverLUsparse(H,g,A,b);
    elseif strcmpi(solver,"LDLdense")
        [x,lambda] = EqualityQPSolverLDLdense(H,g,A,b);
    elseif strcmpi(solver,"LDLsparse")
        [x,lambda] = EqualityQPSolverLDLsparse(H,g,A,b);
    elseif strcmpi(solver,"RangeSpace")
        %H = full(H); g=full(g); A=full(A); b=full(b);
        [x,lambda] = EqualityQPSolverRangeSpace(H,g,A,b);
    elseif strcmpi(solver,"NullSpace")
        %H = full(H); g=full(g); A=full(A); b=full(b);
        [x,lambda] = EqualityQPSolverNullSpace(H,g,A,b);
    end
end

function [x,lambda] = EqualityQPSolverLUdense(H,g,A,b);
    S = full([H ,-A;-A',sparse(size(b,1),size(b,1))]);
    F = full(-[g;b]);
    [L,U,P] = lu(S);
    xL = U\(L\(P*F));
    [x,lambda] = deal(xL(1:length(g)),xL(length(g)+(1:length(b))));
end

function [x,lambda] = EqualityQPSolverLUsparse(H,g,A,b);
    S = [H ,-A;-A',sparse(size(b,1),size(b,1))];
    F = -[g;b];
    [L,U,P] = lu(S);
    xL = U\(L\(P*F));
    [x,lambda] = deal(xL(1:length(g)),xL(length(g)+(1:length(b))));
end

function [x,lambda] = EqualityQPSolverLDLdense(H,g,A,b);
    S = full([H ,-A;-A',sparse(size(b,1),size(b,1))]);
    F = full(-[g;b]);
    [L,D,P] = ldl(S);
    xL = P*(L'\(D\(L\(P'*F))));
    [x,lambda] = deal(xL(1:length(g)),xL(length(g)+(1:length(b))));
end

function [x,lambda] = EqualityQPSolverLDLsparse(H,g,A,b);
    S = [H ,-A;-A',sparse(size(b,1),size(b,1))];
    F = -[g;b];
    [L,D,P] = ldl(S);
    xL = P*(L'\(D\(L\(P'*F))));
    [x,lambda] = deal(xL(1:length(g)),xL(length(g)+(1:length(b))));
end


function [x,lambda] = EqualityQPSolverRangeSpace(H,g,A,b);
    L = chol(H);
    v = L\(L'\g);
    aHa = A'*(L\(L'\A));
    lambda = aHa\(b+A'*v);
    x = L\(L'\(A*lambda-g));
end

function [x,lambda] = EqualityQPSolverNullSpace(H,g,A,b);
    n = size(A,1);
    ma = size(A,2);
    [Q,R] = qr(A);
    R = R(1:ma,:);
    Q1 = Q(:,1:ma);
    Q2 = Q(:,(ma+1):n);
    xy = R'\b;
    t1 = (Q2'*H*Q2);
    t2 = -Q2'*(H*Q1*xy+g);
    xz = t1\t2;
    x = Q1*xy+Q2*xz;
    lambda = R\(Q1'*(H*x+g));
end
