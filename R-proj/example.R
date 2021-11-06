library(volesti)

A = R.matlab::readMat('A.mat')[[1]]
b = c(R.matlab::readMat('b.mat')[[1]])
V = R.matlab::readMat('V.mat')[[1]]
x = c(R.matlab::readMat('x.mat')[[1]])
x0 = c(R.matlab::readMat('x0.mat')[[1]])
Vnorms = c(R.matlab::readMat('Vnorms.mat')[[1]])

N=1000
walk_length = 5

b2 = b - A%*%x0
samples = sample_component(A, b2, x-x0, N, walk_length, x0, V, Vnorms, 1, 1)

