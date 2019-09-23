library(volesti)
library(R.matlab)

modelmat = readMat('model.mat')
modelmat = modelmat$model

Aeq = modelmat[1]
Aeq = Aeq[[1]]
Aeq=matrix(Aeq,ncol = ncol(Aeq), nrow = nrow(Aeq))

lb = modelmat[6]
lb = lb[[1]]

ub = modelmat[7]
ub = ub[[1]]

beq = modelmat[3]
beq = beq[[1]]

c = modelmat[8]
c = c[[1]]

d = dim(Aeq)[2]
d = length(lb)

A = rbind(diag(d), -diag(d))
b = rbind(ub, -lb)


points = SamplePoints(A = A, b = b, Aeq = Aeq, beq = beq, random_walk = "CDHR",walk_length = 5)
