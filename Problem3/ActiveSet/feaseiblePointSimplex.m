function x0 = feaseiblePointSimplex(A, b)

    [m, n] = size(A);


    %make x0 = [x; t; s1; s2]
    x = zeros(n, 1);
    t = max(abs(b));
    s1 = t*ones(m, 1) - b;
    s2 = t*ones(m, 1) + b;
    x0 = [x; t; s1; s2];

    %create g so that it minimizes t only
    g = [zeros(n, 1); 1; zeros(2*m, 1)];

    %create A to satisfy the new equality conditions
    % Ax + t >= +-b     Ax + t - s = +-b
    A = [A, ones(m, 1), -eye(m), zeros(m);
        -A, ones(m, 1), zeros(m), -eye(m)];

    %find the feasible point using simplex
    x0 = Simplex(g, A, x0);
    t = x0(n+1);
    x0 = x0(1:n);

    if t > 1e-6
        error("Infeasible")
    end
end
