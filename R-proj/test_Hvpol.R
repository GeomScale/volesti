library(volesti)

#P=GenCross(10,'V')
d=40
k=80
P=GenRandVpoly(d,k, body = "cube")

vol = volume(P, algo = "CB", random_walk = "BilW", parameters = list("hpoly"=TRUE, "nfacets"=10*d))

print(vol)

#vol = volume(P, algo = "CB",rounding = TRUE, random_walk = "BilW")

#print(vol)

#vol = volume(P, algo = "CB", rounding = TRUE, random_walk = "BilW")

#print(vol)
