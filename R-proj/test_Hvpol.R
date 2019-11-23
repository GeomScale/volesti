library(volesti)

#P=GenCross(10,'V')
d=80
k=160
P=GenRandVpoly(d,k, body = "cube")

vol = volume(P, algo = "CB", random_walk = "BilW", parameters = list("hpoly"=TRUE, "nfacets"=5*ceiling(log2(d))*d))

print(vol)

#vol = volume(P, algo = "CB",rounding = TRUE, random_walk = "BilW")

#print(vol)

#vol = volume(P, algo = "CB", rounding = TRUE, random_walk = "BilW")

#print(vol)
