cvx_begin
    variable a(2)
    variable b
    minimize (norm(a) + gamma*sum(pos(1 - X*a + b) + pos(1 + Y*a - b)))
cvx_end

x = [a;b];