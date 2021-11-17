library(volesti)

A = R.matlab::readMat('A.mat')
A = A[[1]]

b = R.matlab::readMat('b.mat')
b = b[[1]]

x = R.matlab::readMat('x.mat')
x = x[[1]]

x0 = R.matlab::readMat('x0.mat')
x0 = x0[[1]]

X = R.matlab::readMat('X.mat')
X = X[[1]]

V = R.matlab::readMat('V.mat')
V = V[[1]]

Sind = R.matlab::readMat('Sind.mat')
Sind = Sind[[1]]


res = get_max_cap(X, A, b, x0, V, Sind, 2000)

mu = res$center
sigma = cov(t(X))

Y0 = sample_gaussian(A,b,mu,mu,sigma,0.001,20,1000,5,x0)

R.matlab::writeMat(Y0=Y0,'Y0.mat')

vol = compute_component_volume(A,b,mu,sigma,x0,X,10,0.05)
print(vol)


