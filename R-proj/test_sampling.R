library(volesti)

P = gen_rand_hpoly(2, 7)
p = sample_points(P, n=2000, distribution = list("density"="gaussian", "variance"=1), random_walk = list("walk"="ExactHMC"))
