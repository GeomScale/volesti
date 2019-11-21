library(volesti)

#P=GenCross(10,'V')
P=GenRandVpoly(30,60)

vol = volume(P, algo = "CB", random_walk = "BilW",rounding = TRUE, parameters = list("hpoly"=TRUE, "nfacets"=120))

print(vol)

vol = volume(P, algo = "CB", random_walk = "BilW", rounding = TRUE)

print(vol)
