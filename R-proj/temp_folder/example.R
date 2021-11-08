library(volesti)

A = R.matlab::readMat('A.mat')[[1]]
b = c(R.matlab::readMat('b.mat')[[1]])
V = R.matlab::readMat('V.mat')[[1]]
x = c(R.matlab::readMat('x.mat')[[1]])
x0 = c(R.matlab::readMat('x0.mat')[[1]])
Vnorms = c(R.matlab::readMat('Vnorms.mat')[[1]])
sigma = R.matlab::readMat('sigma.mat')[[1]]
c = R.matlab::readMat('c.mat')[[1]]
c=c[1]

N=20000
walk_length = 1

b2 = b - A%*%x0
samples = sample_component(A, b, x, N, walk_length, x0, V)
q=psrf_univariate(samples)
print(q)

X = estimate_component(A, b, x, 1, 2000, x0,V,1.01,1,0.1,1,1200,TRUE)

Y = check_convergence(A,b,V,x0,samples,10,0.1,0.1,FALSE)

print(Y)

y = return_first_inside(A, b, V, x0, samples)
