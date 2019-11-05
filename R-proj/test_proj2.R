library(volesti)

d = 3
k = 7
m = 20
P = GenRandZonotope(d, k, dist = "uniform")

T = P$G
P2 = GenRandHpoly(k,m)
A = P2$A
b = P2$b

P3 = IntPoly$new(T=t(T), A=A, b=b)

#vol1 = volume(P, algo = "CB", random_walk = "BilW")

vol2 = volume(P3, algo = "CB", random_walk = "RDHR")
