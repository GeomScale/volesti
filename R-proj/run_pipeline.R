library(volesti)

sigma = R.matlab::readMat('./data/sigma_i_1.mat')[[1]]
c = R.matlab::readMat('./data/c_i_1.mat')[[1]]
c=c[1]

d = dim(sigma)[2]-1

parameters = get_parameters(d)
M = 10000

samples = sample_ptfs_constant_volatility(sigma, c, M, parameters)

correctness = check_correctness(samples, sigma, c)
print(correctness)

