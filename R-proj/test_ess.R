library(volesti)

P = gen_rand_hpoly(4,13)
p = sample_points(P,n=24)
ess(p)