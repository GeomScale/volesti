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

mu = R.matlab::readMat('mu.mat')
mu = mu[[1]]

d = 3

counter = 0
N=100000
Y = boundary_randsphere(d, N) + kronecker(matrix(1, 1, N), matrix(x0, ncol = 1))

for (i in 1:N) {
  y = Y[,i]
  if(is_in_component(y, x0, A, b, V, Sind, TRUE)){
    counter = counter + 1
  }
}
log_vol = log_volume_n_sphere(d) + log(counter/N)
vol = exp(log_vol)
#res = get_max_cap(X, A, b, x0, V, Sind, 2000)

#mu = res$center
sigma = cov(t(X))

Y0 = sample_fischer(A,b,mu,mu,50,20,1000,5,x0)

R.matlab::writeMat(Y0=Y0,'Y0.mat')

vol = compute_fischer_volume(A,b,mu,sigma,x0,X,10,0.05)
log_int = log_integral_fischer(d, vol$last_variance)
print(vol$ratio*exp(log_int))

res1 = compute_fischer_annealing(A,b,mu,sigma,x0,X,10)
print(res1)

res2 = compute_fischer_ratios(A,b,mu,sigma,x0,X,res1$variances, res1$ratios, 1004, 10, 0.1)
print(res2)
print(exp(res2))

