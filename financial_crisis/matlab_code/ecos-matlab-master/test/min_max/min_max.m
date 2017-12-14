cvx_begin
    variable x(3)
    maximize (min(x(1) - x(2), x(3)))
    subject to
        x(2) == a*x(1)
        x(3) == a*x(2)
        abs(x(2)) <= b
        max(x(1),x(3)) <= b
cvx_end