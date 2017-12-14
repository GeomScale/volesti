c_ = sparse(8,1);
c_(1) = 1;
b_ = sparse(7,1);
b_(3:3) = b;
A_ = sparse(7, 8);
A_(1:1, 4:4) = 1.0; A_(1:1, 5:5) = -1.0;
A_(2:2, 6:7) = A; A_(2:2, 4:4) = -1.0;
A_(3:3, 5:5) = 1.0;
A_(4:5, 6:7) = 1.0*speye(2, 2); A_(4:5, 2:3) = -1.0*speye(2, 2); A_(4:5, 8:8) = -1.0*ones(2, 1);
A_(6:6, 8:8) = 1.0;
A_(7:7, 6:7) = ct; A_(7:7, 1:1) = -1.0;
cvx_begin
variable x_codegen(8)
minimize (c_'*x_codegen)
A_*x_codegen == b_
x_codegen(2:3) >= 0
cvx_end
t0 = x_codegen(1:1);
t3 = x_codegen(2:3);
t1 = x_codegen(4:4);
pb = x_codegen(5:5);
x = x_codegen(6:7);
c0d0 = x_codegen(8:8);
ecos_optval = 1*cvx_optval;
