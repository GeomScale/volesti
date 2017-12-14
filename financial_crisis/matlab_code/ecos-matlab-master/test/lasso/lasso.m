cvx_begin
    variable x(n)
    minimize (sum_square(A*x - b) + lambda*norm(x,1))
    subject to
        -0.1 <= x <= 0.1
cvx_end