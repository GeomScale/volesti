% very simple MPC problem
 
% solve with cvx
cvx_begin
    variables x0 u0 x1 u1 x2 u2 x3
    minimize (q*square(x0) + r*square(u0) + q*max(x1,1) - r*min(u1,0) + quad_over_lin(x2,1/q) - r*sqrt(u2) + p*inv_pos(x3))
    subject to 
       x1 == a*x0 + b*u0;
       x2 == a*x1 + b*u1;
       x3 == a*x2 + b*u2;
       x0 == x_0;
       0 <= u0 <= 1;
       0 <= u1 <= 1;
       0 <= u2 <= 1;
cvx_end

x = [x0;x1;x2;x3;u0;u1;u2];