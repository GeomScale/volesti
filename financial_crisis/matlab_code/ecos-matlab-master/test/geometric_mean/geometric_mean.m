cvx_begin
    variable x(2)
    maximize (x(1) + geo_mean(x))
    subject to
        a <= x(1);
        x(1) <= b
        x(2) <= x(1)
cvx_end