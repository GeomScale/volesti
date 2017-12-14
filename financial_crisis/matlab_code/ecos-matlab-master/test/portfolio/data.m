n = 10;
m = 3;

F = randn(n,m);
Ft = F';
d = sqrt(rand(n,1));
D = diag(d);

B = 20;

mu = rand(n,1);
mut = mu';

gamma = 1;

params.F = F; params.d = d; params.B = B; params.mu = mu; params.gamma = gamma;
params.D = D;
