d = length(x0)
X = boundary_randsphere(d, N) + kronecker(matrix(1, 1, N), matrix(x0, ncol = 1))

for (i in 1:N) {
  x=X[,i]
  if(is_in_component(x,x0,A,b,V,S_Vindices[[1]], TRUE)){
    break
  }
}

Y = sample_component_psrf(A,b,x,10000,5,x0,1.01)
mu = rowMeans(Y)
mu = mu / sqrt(sum(mu^2))
sigma = cov(t(Y))

samples = sample_gaussian(A,b,x,mu,sigma,1,5,1000,5,x0)

Sind = S_Vindices[[1]]

R.matlab::writeMat(X=samples,'samples.mat')
R.matlab::writeMat(A=A,'A.mat')
R.matlab::writeMat(b=b,'b.mat')
R.matlab::writeMat(mu=mu,'mu.mat')
R.matlab::writeMat(sigma=sigma,'sigma.mat')
R.matlab::writeMat(x=x,'x.mat')
R.matlab::writeMat(x0=x0,'x0.mat')
R.matlab::writeMat(V=V,'V.mat')
R.matlab::writeMat(Sind=Sind,'Sind.mat')

