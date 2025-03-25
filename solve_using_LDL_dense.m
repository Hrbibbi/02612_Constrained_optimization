function [x]=solve_using_LDL(A,b)
[L,D,P] = ldl(A);
x = P*(L'\(D\(L\(P'*b))));

end