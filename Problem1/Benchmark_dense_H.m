clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem 1: Equality Constrained Convex QP (subproblem 4-5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rng(42069); % Set seed for reproducability while debugging

npoint = 10;
ProblemSizes = round(logspace(2,3,npoint));  % Problem sizes
beta = 0.8;                     % Constraint to design variable density
alpha = 0.1;                     % Shift value
nrep = 10;                        % Number of repitions

% Time containers
tLUdens = zeros(npoint,1);
tLUsparse = zeros(npoint,1);
tLDLdens = zeros(npoint,1);
tLDLsparse = zeros(npoint,1);
tRangeSpace = zeros(npoint,1);
tNullSpace = zeros(npoint,1);

for j= 1:nrep
    disp(j)
    for i=1:npoint
        n = ProblemSizes(i);  % Numer of optimization variables
        m = round(beta*n);    % Number of Constraints
    
        % Generate test problem:
        Aeq = sprandn(n,m,0.15);
        M = sprandn(n,n,0.15);
        H = M*M'+alpha*speye(n);
        x =  randn(n,1);
        lambda = randn(m,1);
        beq = Aeq'*x;
        g = randn(n,1);
    
        % Matlab Build in solver:
        %xM = quadprog(H,g,[],[],Aeq',beq);
    
        tic;
        [x1,lambda1] = EqualityQPSolver(H,g,Aeq,beq,'LUdense');
        tLUdens(i) =  tLUdens(i)+toc;
        tic; 
        [x2,lambda2] = EqualityQPSolver(H,g,Aeq,beq,'LUsparse');
        tLUsparse(i) = tLUsparse(i)+ toc;
        tic;
        [x3,lambda3] = EqualityQPSolver(H,g,Aeq,beq,'LDLdense');
        tLDLdens(i) = tLDLdens(i)+toc;
        tic; 
        [x4,lambda4] = EqualityQPSolver(H,g,Aeq,beq,'LDLsparse');
        tLDLsparse(i) = tLDLsparse(i)+toc;
        tic; 
        [x5,lambda5] = EqualityQPSolver(H,g,Aeq,beq,'RangeSpace');
        tRangeSpace(i) = tRangeSpace(i)+ toc;
        tic; 
        [x6,lambda6] = EqualityQPSolver(H,g,Aeq,beq,'NullSpace');
        tNullSpace(i) = tNullSpace(i)+ toc;
    end
end
figure()
loglog(ProblemSizes,tLUdens/nrep,'--o','DisplayName','LU dense')
hold on;
loglog(ProblemSizes,tLUsparse/nrep,'--^','DisplayName','LU sparse')
loglog(ProblemSizes,tLDLdens/nrep,'--v','DisplayName','LDL dense')
loglog(ProblemSizes,tLDLsparse/nrep,'-->','DisplayName','LDL sparse')
loglog(ProblemSizes,tRangeSpace/nrep,'--<','DisplayName','Range Space')
loglog(ProblemSizes,tNullSpace/nrep,'--square','DisplayName','Null Space')
%title(['beta = ' num2str(beta)])
legend('Location','northwest','Interpreter','latex')
xlabel("Problem size $n$ [-]",'Interpreter','latex')
ylabel("Computation time [s]",'Interpreter','latex')
exportgraphics(gcf,['Problem1_benchmark_dense_beta_' num2str(beta) '.png'],'Resolution',300)

