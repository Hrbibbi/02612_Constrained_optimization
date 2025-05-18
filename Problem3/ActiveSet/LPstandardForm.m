function [g_std, A_std, b_std] = LPstandardForm(g, A, b, u)
    %convert the lp problem into standard form

    [m, n] = size(A);

    %now we update the variables to be in standard form
    g_std = [g; zeros(m, 1)];

    A_std = [A', zeros(n, m)
             eye(m), eye(m)];

    b_std = [b;u];
end

