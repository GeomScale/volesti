n = 2;
A = randn(1,n);
b = rand();
c = [1;1];

% a1 = A(1,1);
% a2 = A(1,2);
% c1 = c(1);
% c2 = c(2);
% zero = 0;

ct = c';

params.A = A;
params.b = b;
params.c = c;
params.ct = ct;