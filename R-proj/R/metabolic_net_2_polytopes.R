metabolic_net_2_polytope <- function(path) {

modelmat = R.matlab::readMat(path)
modelmat = modelmat[1]
modelmat = modelmat[[1]]

el = 11
if(dim(modelmat[[11]])[1] > dim(modelmat[[11]])[2]) {
  el = el -1
}

Aeq = modelmat[[el]]

#Aeq=matrix(Aeq,ncol = ncol(Aeq), nrow = nrow(Aeq))

lb = as.vector(modelmat[[el+1]])
ub = as.vector(modelmat[[el+2]])
beq = as.vector(modelmat[[el+3]])
c = as.vector(modelmat[[el+4]])

d = dim(Aeq)[2]

A = rbind(diag(d), -diag(d))
b = c(ub, -lb)

HP = Hpolytope$new(A = A, b = b, Aeq = Aeq, beq = beq)
return(HP)
}

