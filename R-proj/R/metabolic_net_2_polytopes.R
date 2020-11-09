metabolic_net_2_polytope <- function(path) {

modelmat = R.matlab::readMat(path)
modelmat = modelmat[1]
modelmat = modelmat[[1]]

Aeq = modelmat[[11]]

#Aeq=matrix(Aeq,ncol = ncol(Aeq), nrow = nrow(Aeq))

lb = as.vector(modelmat[[12]])
ub = as.vector(modelmat[[13]])
beq = as.vector(modelmat[[14]])
c = as.vector(modelmat[[15]])

d = dim(Aeq)[2]

A = rbind(diag(d), -diag(d))
b = c(ub, -lb)

HP = Hpolytope$new(A = A, b = b, Aeq = Aeq, beq = beq)
return(HP)
}

