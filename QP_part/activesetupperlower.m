function [x, lambda_l, lambda_u, lambda_Cneg, lambda_C] = activesetupperlower(H,g,C,l,u,dl,du)
    %---------------------------
    %     Initalization step
    %---------------------------
    n = length(l);
    p = size(C,2);
    A = [ eye(n),      -eye(n),    C,     -C     ];
    
    %b set to be consistent with the Ax+b=0 form 
    %instead of the usual Ax=b form
    b = [ -l;           u;        -dl;     du   ];
    
    n_constraints = size(A,2);
    lambda_full = zeros(n_constraints, 1);

    % start at lower‐bound corner:
    x = l;
    Workingset = 1:n;
    tol = 1e-4;
    %---------------------------------
    %           Iteration step
    %---------------------------------
    for k = 1:1000
        %---------------------------------
        %           KKT system solve
        %---------------------------------
        Ak   = A(:,Workingset);
        gk   = H*x + g;
        KKT  = [ H,    -Ak; -Ak', zeros(length(Workingset)) ];
        rhs  = [ -gk; zeros(length(Workingset),1) ];
        sol  = KKT \ rhs;
        p    = sol(1:n);
        lambda = sol(n+1:end);
        
        %---------------------------------
        %       step check
        %---------------------------------
        if norm(p) < tol
            if all(lambda >= 0)
            %Converged
            lambda_full(:) = 0;
            lambda_full(Workingset) = lambda;
            break
            else
            %Constraint update
            [~, index] = min(lambda);
            Workingset(index)   = [];
            end
        else
            %--------------------------------------
            %           alpha computation
            %-------------------------------------- 
            cinactive   = setdiff(1:n_constraints, Workingset);
            Ap     = A(:,cinactive)' * p;
            mask   = (Ap < 0);
            cinactive = cinactive(mask);
            alphas = -( A(:,cinactive)'*x + b(cinactive) ) ./ Ap(mask);
            [alpha, idx] = min(alphas);
            if alpha < 1
                x     = x + alpha*p;
                Workingset     = [Workingset, cinactive(idx)];
            else
                x = x + p;
            end
        end
    end
    lambda_l     = lambda_full(1:200);            % x ≥ l
    lambda_u     = lambda_full(200+1:400);        % x ≤ u
    lambda_C     = lambda_full(401:600);    % C^T x ≤ du
    lambda_Cneg  = lambda_full(601:end);    % C^T x ≥ dl
end
