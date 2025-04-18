% Initial guess and active set
x = [4;1
    ]; wactive = [];
n_constraint = size(b,1);

for k = 0:10
    fprintf('iteration: %3i , x1 = %7.3f , x2 = %7.3f, active constraints:', k, x(1),x(2))
    if isempty(wactive)
        fprintf('\n')
    else
        disp(wactive)
    end
    [p,lambda] = SolveQP(H,g,A,b,x,wactive); %Solve KKT
    if all(abs(p)<1e-10) % Check if converged or remove constraint
        if all(lambda>=0)
            break
        else
            wactive(find(lambda == min(lambda)))=[];
        end
    else
        winactive = setdiff(1:n_constraint,wactive);
        winactive = winactive(A(winactive,:)*p<0);
        alphas = -(A(winactive,:)*x+b(winactive))./(A(winactive,:)*p);
        alpha = min(alphas);
        if alpha<1
            xnew = x+alpha*p;
            wactive = union(wactive,winactive(alpha==min(alphas))); % Add constraint to active set
        else
            xnew = x+p;
        end
        plot([x(1),xnew(1)],[x(2),xnew(2)],'-or')
        x = xnew;
    end
end