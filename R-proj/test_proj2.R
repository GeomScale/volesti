library(volesti)

d =10
k = 20
m = 40
P = GenRandZonotope(d, k, dist = "uniform")

T = P$G
P2 = GenRandHpoly(k,m)
ball = InnerBall(P2)
print(ball)
A = P2$A
b = P2$b

P3 = IntPoly$new(T=t(T), A=A, b=b)

#vol1 = volume(P, algo = "CB", random_walk = "BilW")

#vol1 = volume(P3, algo = "CB", random_walk = "RDHR")

vol2 = volume(P3, algo = "CB", random_walk = "BilW", rounding= FALSE)

#print(vol1)

print(vol2)
