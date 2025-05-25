function [x] = IpoptWrapper(x0in);

global iter_data;
iter_data = [x0in'];

% Set the IPOPT options
options.ipopt.tol = 1e-12;
options.ipopt.hessian_approximation = 'exact';
options.ipopt.derivative_test = 'second-order';
options.ipopt.print_level = 5;
% Number of variables and constraints
n = 2;
m = 2;

% Initial guess
x0 = x0in;

% Upper and lower bounds on optimization variables and constraints
options.lb = [-4;-2.5];
options.ub = [1;4];
options.cl = [0,18];
options.cu = [0,34];

% Initialize the dual point
options.zl = ones(n, 1);
options.zu = ones(n, 1);
options.lambda = ones(m, 1);

% The callback functions
funcs.objective = @eval_f;
funcs.constraints = @eval_g;
funcs.gradient = @eval_grad_f;
funcs.jacobian = @eval_jac_g;
funcs.jacobianstructure = @eval_jac_g_sparsity_pattern;
funcs.hessian = @eval_h;
funcs.hessianstructure = @eval_h_sparsity_pattern;

% Run IPOPT
[x, info] = ipopt(x0, funcs, options);

figure(1);
hold on;
plot(iter_data(:,1),iter_data(:,2),':ok');
end

function f = eval_f(x)
    tmp1 = (x(1)^2+x(2)-11);
    tmp2 = x(1)+x(2)^2-7;
    f = tmp1^2+tmp2^2;
    global iter_data;
    iter_data = [iter_data; x'];
end

function grad_f = eval_grad_f(x)
    tmp1 = (x(1)^2+x(2)-11);
    tmp2 = x(1)+x(2)^2-7;
    grad_f(1) = 4*x(1)*tmp1+2*tmp2;
    grad_f(2) = 2*tmp1+4*x(2)*tmp2;
end

function g = eval_g(x)
    % Equality constraints
    tmp = (x(1)+2);
    g(1) = tmp^2-x(2)-3;
    % Inequality constraints 
    g(2) = -4*x(1)+10*x(2);
end



function j = eval_jac_g(x)
    % Row and column indices
    rows = [1, 1, 2, 2];
    cols = [1, 2, 1, 2];
    % Evaluate the Jacobian
    vals = [2*(x(1) + 2), -1, -4, 10];
    % Jacobian
    j = sparse(rows, cols, vals);
end

function j = eval_jac_g_sparsity_pattern()
% Row and column indices
    rows = [1, 1, 2, 2];
    cols = [1, 2, 1, 2];
% Make Jacobian sparse
    j = sparse(rows, cols, 1);
end

function h = eval_h(x, sigmaf, lambda)
% tempary variables
    tmp1 = (x(1)^2+x(2)-11);
    tmp2 = x(1)+x(2)^2-7;
% Row and column indices
    rows = [1, 2, 2];
    cols = [1, 1, 2];
    vals = zeros(3,1);
% Evaluate the hessian of objective
    vals(1) = sigmaf*(4*tmp1+8*x(1)^2+2);
    vals(2) = sigmaf*(4*(x(1)+x(2)));
    vals(3) = sigmaf*(4*tmp2+8*x(2)^2+2);
% and the contribution from the constriants
    vals(1) = vals(1) + lambda(1)*2;
% Hessian
    h = sparse(rows, cols, vals);
end

function h = eval_h_sparsity_pattern()
% Row and column indices
    rows = [1, 2, 2];
    cols = [1, 1, 2];
% Make Jacobian sparse
    h = sparse(rows, cols, 1);
end