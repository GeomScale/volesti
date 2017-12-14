cvx_begin
    variable x(2)
    minimize (-x(1) + x(2) + inv_pos(sqrt(a*x(1) - b*x(2))))
        pos(x(1)) <= a
        pos(-x(2)) <= b
cvx_end