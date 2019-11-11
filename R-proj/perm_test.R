library(volesti)

d =5
k = 20
m = 40
P = Permutaedron$new(3)

vol2 = volume(P, algo = "CB", random_walk = "RDHR")

print(vol2)
