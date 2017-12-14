n = 2;       % dimension
m = 10;    % number of samples
X = randn(m,n) + 1.5;
Y = randn(m,n) - 1.5;

gamma = 1;

o = ones(m,1);

params.X = X;
params.Y = Y;
params.gamma = gamma;
params.o = o;