library(volesti)
library(R.matlab)

modelmat = readMat('roundedPolytope.mat')
roundedPoly = modelmat$poly

A = roundedPoly[1]
A = A[[1]]
A = A[[1]]

b = roundedPoly[2]
b = b[[1]]
b = b[[1]]

T = roundedPoly[3]
T = T[[1]]
T = T[[1]]

p = roundedPoly[4]
p = p[[1]]
p = p[[1]]

G = roundedPoly[5]
G = G[[1]]
G = G[[1]]

shift = roundedPoly[6]
shift = shift[[1]]
shift = shift[[1]]

P = Hpolytope$new(A=A,b=b)
#P2 = RoundPoly(P)$P
points = SamplePoints(P, random_walk = "BilW",walk_length = 100, n = 1000)

#points = T%*%points +  matrix(rep(p,1000),ncol=1000)
points = G%*%points + matrix(rep(shift,1000),ncol=1000)

writeMat(paste0(getwd(),"/samples.mat"),points=points)
