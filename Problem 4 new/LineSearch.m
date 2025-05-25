function [alpha] = LineSearch(xk,deltaX,lambda,mu,obj,ceq,cineq,gl,gu,xl,xu);
% Evaluate objective and constraints
[f,df] = obj(xk);
[h] = ceq(xk);
[g] = cineq(xk);
G = [gu-g;g-gl;xu-xk;xk-xl];

alpha = 1;
c = f+mu'*abs(h)+lambda'*abs(min(0,G));
b = df'*deltaX-mu'*abs(h)-lambda'*abs(min(0,G));

for i = 1:100
    x = xk+alpha*deltaX;

    % Evaluate obj and con
    [f] = obj(x);
    [h] = ceq(x);
    [g] = cineq(x);
    G = [gu-g;g-gl;xu-xk;xk-xl];
    
    phi = f'+mu'*abs(h)+lambda'*abs(min(0,G));

    if phi<=(c+0.1*b*alpha)
        break
    else
        a = (phi-(c+b*alpha))/(alpha^2);
        alphamin = -b/(2*a);
        alpha = min(0.9*alpha,max(alphamin,0.1*alpha));
    end
end