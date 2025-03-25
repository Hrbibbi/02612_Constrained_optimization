function [x,lambda] = solve_using_rangespace(H,g,A,b)
%assumption H is identity
A=A';
v=g';
HA=A'*A;
lambda=HA\(b+A'*v);
x=A*lambda-g';
end