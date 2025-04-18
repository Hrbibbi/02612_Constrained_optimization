function [p,lambda] = SolveQP(H,g,A,b,xk,W)
    Ar = A(:,W);
    br = b(W);
    S = [H ,-Ar;-Ar',zeros(size(br,1),size(br,1))];
    gk = H*xk+g;
    f = [-gk;zeros(size(br))];
    res = S\f;
    p = res(1:size(H,1));
    lambda = res((size(H,1)+1):end);
end