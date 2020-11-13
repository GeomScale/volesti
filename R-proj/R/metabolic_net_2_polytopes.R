metabolic_net_2_polytope <- function(path, remove_biomass = FALSE, Recon2 = FALSE, Recon3 = FALSE) {

modelmat = R.matlab::readMat(path)
modelmat = modelmat[1]
modelmat = modelmat[[1]]

if (Recon2) {
  
  Aeq = modelmat[[1]]
  
  #Aeq=matrix(Aeq,ncol = ncol(Aeq), nrow = nrow(Aeq))
  
  lb = as.vector(modelmat[[3]])
  ub = as.vector(modelmat[[4]])
  beq = as.vector(modelmat[[21]])
  obj = as.vector(modelmat[[6]])
  
  if (remove_biomass) {
    biomass_col = which(obj == 1)
    Aeq = Aeq[, -c(biomass_col)]
    A = A[, -c(biomass_col)]
    ub = ub[-c(biomass_col)]
    lb = lb[-c(biomass_col)]
  }
  
  d = dim(Aeq)[2]
  
  A = rbind(diag(d), -diag(d))
  b = c(ub, -lb)
  
  HP = Hpolytope$new(A = as.matrix(A), b = c(b), Aeq = as.matrix(Aeq), beq = c(beq))
  return(HP)
}

if (Recon3) {
  
  Aeq = modelmat[[1]]
  
  #Aeq=matrix(Aeq,ncol = ncol(Aeq), nrow = nrow(Aeq))
  
  lb = as.vector(modelmat[[6]])
  ub = as.vector(modelmat[[7]])
  beq = as.vector(modelmat[[3]])
  obj = as.vector(modelmat[[8]])
  
  if (remove_biomass) {
    biomass_col = which(obj == 1)
    Aeq = Aeq[, -c(biomass_col)]
    A = A[, -c(biomass_col)]
    ub = ub[-c(biomass_col)]
    lb = lb[-c(biomass_col)]
  }
  
  d = dim(Aeq)[2]
  
  A = rbind(diag(d), -diag(d))
  b = c(ub, -lb)
  
  HP = Hpolytope$new(A = as.matrix(A), b = c(b), Aeq = as.matrix(Aeq), beq = c(beq))
  return(HP)
}

el = 11
if(dim(modelmat[[11]])[1] > dim(modelmat[[11]])[2]) {
  el = el -1
}

Aeq = modelmat[[el]]

#Aeq=matrix(Aeq,ncol = ncol(Aeq), nrow = nrow(Aeq))

lb = as.vector(modelmat[[el+1]])
ub = as.vector(modelmat[[el+2]])
beq = as.vector(modelmat[[el+3]])
obj = as.vector(modelmat[[el+4]])

if (remove_biomass) {
  biomass_col = which(obj == 1)
  Aeq = Aeq[, -c(biomass_col)]
  A = A[, -c(biomass_col)]
  ub = ub[-c(biomass_col)]
  lb = lb[-c(biomass_col)]
}

d = dim(Aeq)[2]

A = rbind(diag(d), -diag(d))
b = c(ub, -lb)

HP = Hpolytope$new(A = A, b = b, Aeq = Aeq, beq = beq)
return(HP)

}

