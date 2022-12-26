preprocess_with_quadprog <- function(A, b, Aeq, beq) {
  
  d = dim(A)[2]
  A_store = A

  minWeights = c()
  maxWeights = c()
  
  A = rbind(Aeq,A)
  b= c(beq,b)

  A = -t(A)
  b = -b
  
  Dmat = diag(d) * 1e-10
  
  for (j in 1:d) {
    
    dvec = -A_store[j,]
    ret = quadprog::solve.QP(Dmat, dvec, A, b, meq=dim(Aeq)[1])
    max_dist = dvec %*% ret$solution

    dvec = A_store[j,]
    ret = quadprog::solve.QP(Dmat, dvec, A, b, meq=dim(Aeq)[1])
    min_dist = dvec %*% ret$solution
    
    minWeights = c(minWeights, min_dist)
    maxWeights = c(maxWeights, max_dist)
    
    width = abs(max_dist - min_dist)
    
    if (width < 1e-07) {
      Aeq = rbind(Aeq, A[j,])
      beq = c(beq, min_dist)
    }
  }
  
  ret_list = list()
  ret_list$Aeq = Aeq
  ret_list$beq = beq
  ret_list$minWeights = minWeights
  ret_list$maxWeights = maxWeights
  
  return(ret_list)
}
