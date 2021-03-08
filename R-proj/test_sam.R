library(volesti)

P = gen_rand_hpoly(2, 7)
p = sample_points(P, n=2000, distribution = list("density"="exponential","bias"=c(0.7071068,0.7071068), 
                                                "variance"=1), 
                  random_walk = list("walk"="HMC", "step_size"=1, "number_of_steps"=30))