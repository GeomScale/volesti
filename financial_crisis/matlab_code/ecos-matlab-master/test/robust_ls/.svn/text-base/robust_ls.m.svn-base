cvx_begin
    variable x(n)
    minimize( sum_square( A * x - b ) )
    subject to
        C * x == d
        norm( x, Inf ) <= e 
cvx_end