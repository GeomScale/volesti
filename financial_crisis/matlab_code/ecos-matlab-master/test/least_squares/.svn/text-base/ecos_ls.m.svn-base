randn('state', 0);
rand('state', 0);
data;

cvx_begin
    variable x(n)
    minimize (sum_square(A*x - b))
    subject to
        C*x == d
cvx_end
x_cvx = x;
fstar = cvx_optval;

c_ = sparse(334,1);
c_(334) = 1;
b_ = sparse(317,1);
b_(1:1) = -0.5*ones(1, 1);
b_(2:2) = -0.5*ones(1, 1);
b_(203:302) = b;
b_(313:317) = d;
A_ = sparse(317, 334);
A_(1:1, 334:334) = 0.5*speye(1, 1); A_(1:1, 1:1) = -1*speye(1, 1);
A_(2:2, 334:334) = -0.5*speye(1, 1); A_(2:2, 2:2) = -1*speye(1, 1);
A_(3:102, 104:203) = 1*speye(100, 100); A_(3:102, 204:303) = -1*speye(100, 100); A_(3:102, 4:103) = -1*speye(100, 100);
A_(103:202, 304:323) = A; A_(103:202, 104:203) = -1*speye(100, 100);
A_(203:302, 204:303) = 1.0*speye(100, 100);
A_(303:307, 324:328) = 1*speye(5, 5); A_(303:307, 329:333) = -1*speye(5, 5);
A_(308:312, 304:323) = C; A_(308:312, 324:328) = -1*speye(5, 5);
A_(313:317, 329:333) = 1.0*speye(5, 5);
G_ = sparse(104, 334);
G_(1:1:1, 1:1) = -speye(1, 1);
G_(2:1:2, 2:2) = -speye(1, 1);
G_(3:1:3, 3:3) = -speye(1, 1);
G_(4:1:4, 3:3) = -speye(1, 1);
G_(5:1:104, 4:103) = -speye(100, 100);
h_ = zeros(104, 1);
dims.q = [3,101];
dims.l = 0;
[x_codegen, y_, info_] = conelp(full(c_), G_, h_, dims, A_, full(b_));
t0z0 = x_codegen(1:1);
t0z1 = x_codegen(2:2);
t1 = x_codegen(3:3);
t2 = x_codegen(4:103);
t3 = x_codegen(104:203);
t4 = x_codegen(204:303);
x = x_codegen(304:323);
t5 = x_codegen(324:328);
t6 = x_codegen(329:333);
t0 = x_codegen(334:334);
ecos_optval = 1*info_.pcost;

tt = square(norm(A*x_cvx - b));
[norm([1/2*(1-tt); norm(A*x_cvx - b)]) (1/2*(1+tt))]
[norm([t0z1;t1]) t0z0]      % doesn't satisfy cone constraint ||x|| > t
[t0 cvx_optval square(norm(A*x-b))]    % not equal...

% square(s) <= t
% norm(A*x-b) <= s

info_