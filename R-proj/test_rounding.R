library(volesti)
library(Matrix)
library(Rmosek)

P = gen_rand_hpoly(45,100)
P$b = runif(100, 0, 1)

r = rounding_isotropic(P)

