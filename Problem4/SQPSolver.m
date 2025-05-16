function [xk,it] = SQPSolver(x0,obj,ceq,cineq,cl,cu,xl,xu,bfgs,steptype,zeta);
    converged = false;
    % Max iteration and convergence strictness
    itmax = 40; epsilon = 1e-2;

    % Evaluate objective and constraints:
    xk = x0; xplot=[x0'];
    [f,df,ddf] = obj(x0);
    [h,dh,ddh] = ceq(x0);
    [g,dg,ddg] = cineq(x0);
    
    % Number of variables and constraints
    n = length(x0);
    mE = length(h);
    mI = length(g);
    
    % Setup lagrange multipliers:
    mu_k = zeros(mE,1); %Equality
    lambda_k = zeros(2*n+2*mI,1); %Inequality and box
    
    % Compute Hessian
    H = ddf;
    % Add equality constraints
    for i = 1:mE
        H = H-mu_k(i)*ddh(:,:,i);
    end
    % Add inequality constraint
    for i = 1:mI
        H = H-lambda_k(i)*ddg(:,:,i);
        H = H+lambda_k(i+mI)*ddg(:,:,i);
    end

    % BFGS constraint check
    if bfgs == true
        Bk = H;
        Bk = diag([1,1]);
    end
    
    
    for it = 1:itmax
        %% Setup inequality QP
        if bfgs == true
            H = Bk;
        end
        
        dG = [dg';-dg';eye(n);-eye(n)];
        G = [g-cu;-g+cl;xk-xu;-xk+xl];
    
     %% Solve inequality QP
        [pk,mu,lambda] = InteriorPoint(H,df,-dh,h,-dG',G,xk);
        [pk1,~,~,~,tmp] = quadprog(H,df,dG,-G,dh',-h,[],[],xk);
        mu1 = tmp.eqlin;
        lambda1 = tmp.ineqlin;
        pmu = mu-mu_k;
        plambda = lambda-lambda_k;
        
     %% Linesearch
        if strcmpi(steptype,'LineSearch')
            [alpha] = LineSearch(xk,pk,lambda,mu,obj,ceq,cineq,cl,cu,xl,xu);
        elseif strcmpi(steptype,'TrustRegion')
            alpha = min(zeta/norm(pk,inf),1);
        end
     %% Perform step
        lambda_k = lambda_k+alpha*plambda;
        mu_k = mu_k+alpha*pmu;
        
        if bfgs == true
            DL1 = df;
            % Equality Constraints
            for i = 1:mE
                DL1 = DL1+mu_k(i)*dh(:,:,i);
            end
            % Inequality constraint
            for i = 1:mI
                DL1 = DL1-lambda_k(i)*dg(:,:,i);
                DL1 = DL1+lambda_k(i+mI)*dg(:,:,i);
            end
            % box constriaints
            for i = 1:n
                vec = zeros(n,1); vec(i) = 1;
                DL1 = DL1-lambda_k(i+2*mI)*vec;
                DL1 = DL1+lambda_k(i+2*mI+n)*vec;
            end
        end
    
        xk = xk+alpha*pk;
    
    %% Update gradients and hessian
        if bfgs == true
            [f,df] = obj(xk);
            [h,dh] = ceq(xk);
            [g,dg] = cineq(xk);
        else
            [f,df,ddf] = obj(xk);
            [h,dh,ddh] = ceq(xk);
            [g,dg,ddg] = cineq(xk);
        end
        % Lagrangian
        DL2 = df;
        % Equality Constraints
        for i = 1:mE
            DL2 = DL2+mu_k(i)*dh(:,:,i);
        end
        % Inequality constraint
        for i = 1:mI
            DL2 = DL2-lambda_k(i)*dg(:,:,i);
            DL2 = DL2+lambda_k(i+mI)*dg(:,:,i);
        end
        % box constriaints
        for i = 1:n
            vec = zeros(n,1); vec(i) = 1;
            DL2 = DL2-lambda_k(i+2*mI)*vec;
            DL2 = DL2+lambda_k(i+2*mI+n)*vec;
        end 

        if bfgs == true
            qk = DL2-DL1;
 
            %Update BFGS
            temp = (Bk*pk);
            temp2 = pk'*temp;
            if (pk'*qk >= (0.2*temp2))
                theta = 1;
            else
                theta = 0.8*temp2/(temp2-pk'*qk);
            end
            rk = theta*qk+(1-theta)*temp;
            Bk = Bk + rk*(rk')/(pk'*rk)-(temp*(temp'))/(pk'*temp);
        else
            % Compute Hessian
            H = ddf;
            % Add equality constraints
            for i = 1:mE
                H = H-mu_k(i)*ddh(:,:,i);
            end
            % Add inequality constraint
            for i = 1:mI
                H = H-lambda_k(i)*ddg(:,:,i);
                H = H+lambda_k(i+mI)*ddg(:,:,i);
            end
        end
        xplot = [xplot;xk'];
        %% Check KKT
        if (norm(DL2,inf) <= epsilon && norm(h,inf)<=epsilon && all([g-cu;-g+cl;xk-xu;-xk+xl]<=epsilon))
            converged = true;
            break
        end
        figure(1)
        hold on;
        plot(xplot(:,1),xplot(:,2),':ok')
    end
    if converged == true
        fprintf("converged in %2i \n",it)
    else
        fprintf("did not converge \n")
end
