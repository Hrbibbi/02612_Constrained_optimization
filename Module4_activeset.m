clear all; close all;

function [p,lambda] = SolveQP(H,g,A,b,xk,W)
    Ar = A(W,:);
    br = b(W);
    S = [H ,-Ar';-Ar,zeros(size(br,1),size(br,1))];
    gk = H*xk+g;
    f = [-gk;zeros(size(br))];
    res = S\f;
    p = res(1:size(H,1));
    lambda = res((size(H,1)+1):end);
end

function contourplt(H,g,A,b,lim,npoint);
    [x1,x2] = meshgrid(linspace(lim(1),lim(2),npoint),linspace(lim(3),lim(4),npoint));
    x = [x1(:),x2(:)]';
    for i =1:size(x,2)
        fval(:,i) = g'*x(:,i)+0.5*x(:,i)'*H'*x(:,i);
    end
    fval = reshape(fval,size(x1));
    x1 = reshape(x(1,:),size(x1));
    x2 = reshape(x(2,:),size(x1));
    figure()
    contour(x1,x2,fval)
    hold on;
    for i = 1:length(b)
        if A(i,2) ==0
            x2 = linspace(lim(3),lim(4),npoint);
            x1 = -(A(i,2)*x2+b(i))/A(i,1);
            plot(x1,x2,':k');
        else
            x1 = linspace(lim(1),lim(2),npoint);
            x2 = -(A(i,1)*x1+b(i))/A(i,2);
            plot(x1,x2,':k');
        end
    end
    xlim([lim(1:2)])
    ylim([lim(3:4)])
end

% QP problem
H = 2*eye(2);
g = [-2;-5];
A = [1 -2; -1 -2; -1 2; 1 0; 0 1];
b = [2;6;2;0;0];

contourplt(H,g,A,b,[0,4,0,4],20)

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