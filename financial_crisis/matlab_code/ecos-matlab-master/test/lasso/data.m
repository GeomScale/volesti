m = 2; n = 3;
A = randn(m,n);
b = randn(m,1);
lambda = 0.5;
u = 0.1;
l = -0.1;

a11 = A(1,1);
a12 = A(1,2);
a13 = A(1,3);
a21 = A(2,1);
a22 = A(2,2);
a23 = A(2,3);
b1 = b(1);
b2 = b(2);

params.a11 = a11;
params.a12 = a12;
params.a13 = a13;
params.a21 = a21;
params.a22 = a22;
params.a23 = a23;
params.b1 = b1;
params.b2 = b2;

params.u = u;
params.l = l;

params.lambda = lambda;

params.A = A;
params.b = b;