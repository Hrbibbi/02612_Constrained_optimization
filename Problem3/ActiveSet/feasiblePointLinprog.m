function x0 = feasiblePointLinprog(A, b)
    [m, n] = size(A);
    f = [zeros(n,1); 1];
    A_lp = [ A, -ones(m,1);
            -A, -ones(m,1) ];
    b_lp = [ b; -b ];

    lb = [zeros(n,1); 0];

    z = linprog(f, A_lp, b_lp, [], [], lb, []);
    x0 = z(1:n);
    t_opt = z(end);

    if t_opt > 0
        error("not feasible")
    end
end