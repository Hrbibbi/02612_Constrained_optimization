function [x,lambda] = solve_using_nullspace(H,g,A,b)
[Q,Rbar]=qr(A);
n=size(A,1);
m1=size(Rbar,2);
Q1= Q(:,1:m1);
Q2= Q(:,m+1:n);
R=Rbar(1:m1,1:m1);

X_y=R\b;
intermediate=Q2'*H*Q2;
intermediate_2=-Q2'*(H*Q1*X_y+g);
X_z=intermediate\intermediate_2;
x=Q1*X_y+Q2*X_z;
lambda=R\(Q1'*(H*x+g));
end