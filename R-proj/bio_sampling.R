library(volesti)
library(R.matlab)

modelmat = readMat('e_coli_core.mat')
modelmat = modelmat$e.coli.core

Aeq = modelmat[12]
Aeq = Aeq[[1]]

lb = modelmat[13]
lb = lb[[1]]

ub = modelmat[14]
ub = ub[[1]]

beq = modelmat[15]
beq = beq[[1]]

c = modelmat[16]
c = c[[1]]

d = dim(Aeq)[2]
d = length(lb)

A = rbind(diag(d), -diag(d))
b = rbind(ub, -lb)

points = SamplePoints(A = A, b = b, Aeq = Aeq, beq = beq, random_walk = "RDHR",walk_length = 50)
