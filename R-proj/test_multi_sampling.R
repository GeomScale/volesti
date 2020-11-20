library(volesti)
library(Rmosek)

P = gen_rand_hpoly(30, 100)
max_ball = get_max_inner_ball(P$A, P$b)
A = P$A
b = P$b
N=1000
d = P$dimension
source('~/volume_approximation/R-proj/R/multiphase_sampling.R')
p = multiphase_sampling(A = P$A, b = P$b, 
                                  max_ball = max_ball, n = N,
                                  num_rounding_samples = 10*d, 
                                  max_num_samples = 100*d, rounding = FALSE)

print(dim(p))