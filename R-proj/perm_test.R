library(volesti)

d =4
k = 20
m = 40
P = Permutaedron$new(d)

vol2 = volume(P, algo = "CB", random_walk = "RDHR")

print((d)^(d-2))
print(vol2)
