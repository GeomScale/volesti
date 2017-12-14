c_ = sparse(7,1);
c_(1) = 1;
b_ = sparse(6,1);
b_(5:6) = b;
A_ = sparse(6, 7);
A_(1:2, 4:5) = 1.0*speye(2, 2); A_(1:2, 2:3) = 1.0*speye(2, 2); A_(1:2, 6:7) = -1.0*speye(2, 2);
A_(3:4, 1:1) = A; A_(3:4, 4:5) = -1.0*speye(2, 2);
A_(5:6, 6:7) = 1.0*speye(2, 2);
cvx_begin
variable x_codegen(7)
minimize (c_'*x_codegen)
A_*x_codegen == b_
x_codegen(2:3) >= 0
cvx_end
x = x_codegen(1:1);
t1 = x_codegen(2:3);
t0 = x_codegen(4:5);
pb = x_codegen(6:7);
ecos_optval = 1*cvx_optval;
