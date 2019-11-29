library(volesti)

d = 12
k = 24
m = 70
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
#vol2 = volume(P3, algo = "CB", random_walk = "BilW", rounding= FALSE)
vol2=0
vol3=0
vol4=0
if(ball[k+1]>0){
  tim2 = system.time({ vol2 = volume(P3, algo = "CB", random_walk = "BilW", rounding= FALSE, parameters=list("hpoly"=TRUE,"nfacets"=120)) })
  tim3 = system.time({ vol3 = volume(P3, algo = "CB", random_walk = "BilW", rounding= FALSE) })
  tim4 = system.time({ vol4 = volume(P3, algo = "CB", random_walk = "BilW", rounding= TRUE) })
}
#print(vol1)

print(paste0(vol2[1]," ",vol2[2]," ",vol2[3]," ",vol2[4]," ",vol2[5]," ", as.numeric(tim2[3])))
print(paste0(vol3[1]," ",vol3[2]," ",vol3[3]," ",vol3[4]," ",vol3[5]," ", as.numeric(tim3[3])))
print(paste0(vol4[1]," ",vol4[2]," ",vol4[3]," ",vol4[4]," ",vol4[5]," ", as.numeric(tim4[3])))
