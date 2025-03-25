function [x]=solve_using_LU(M,rhs)
[L,U,P]=lu(M);
y=L\(P*rhs);
x=U\y;
end