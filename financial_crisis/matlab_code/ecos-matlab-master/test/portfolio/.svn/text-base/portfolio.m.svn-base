cvx_begin
    variable x(n)

    maximize (mu'*x - gamma*square_pos(norm(F'*x)) - gamma*square_pos(norm(D'*x)))
    subject to
        sum(x) == B
        x >= 0
cvx_end

