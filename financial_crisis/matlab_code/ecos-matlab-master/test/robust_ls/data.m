randn('state',0);

m = 4; n = 2; p = 1;
A = randn(m,n); b = randn(m,1);
C = randn(p,n); d = randn(p,1); e = 5*rand;

a11 = A(1,1);
a12 = A(1,2);
a21 = A(2,1);
a22 = A(2,2);
a31 = A(3,1);
a32 = A(3,2);
a41 = A(4,1);
a42 = A(4,2);

b1 = b(1);
b2 = b(2);
b3 = b(3);
b4 = b(4);

c11 = C(1,1);
c12 = C(1,2);

params.a11 = a11;
params.a12 = a12;
params.a21 = a21;
params.a22 = a22;
params.a31 = a31;
params.a32 = a32;
params.a41 = a41;
params.a42 = a42;

params.b1 = b1;
params.b2 = b2;
params.b3 = b3;
params.b4 = b4;

params.c11 = c11;
params.c12 = c12;

params.A = A;
params.b = b;
params.C = C;
params.d = d;
params.e = e;