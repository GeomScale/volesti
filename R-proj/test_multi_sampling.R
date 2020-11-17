library(volesti)

P = gen_rand_hpoly(10, 100)

p =  sample_points(P, random_walk = list("walk" = "mBiW", 
                                                "walk_length" = 1), n = 100)

print(dim(p))