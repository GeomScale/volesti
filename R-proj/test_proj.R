library(volesti)

d = 15
k = 30
P = GenRandZonotope(d, k, dist = "uniform")

T = P$G
A = rbind(diag(k), -diag(k))
b = rep(1,2*k)

P2 = IntPoly$new(T=t(T), A=A, b=b)

vol1 = volume(P, algo = "CB", random_walk = "BilW")

vol2 = volume(P2, algo = "CB", random_walk = "RDHR")
