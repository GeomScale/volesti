preprocess_with_mosek <- function(P) {
  
  d = P$dimension
  Aeq = P$Aeq
  beq = P$beq
  
  A = P$A
  b = P$b
  
  m = dim(Aeq)[1]
  
  minFluxes = c()
  maxFluxes = c()
  
  row_ind = c()
  col_ind = c()
  values = c()
  
  for (i in 1:m) {
    for (j in 1:d) {
      if (Aeq[i, j] != 0) {
        row_ind = c(row_ind, i)
        col_ind = c(col_ind, j)
        values = c(values, Aeq[i, j])
      }
    }
  }
  
  prob <- list()
  prob$A <- sparseMatrix(row_ind, col_ind, x=values)
  
  # Bound values for constraints
  prob$bc <- rbind(blc=beq, 
                   buc=beq)
  
  # Bound values for variables
  prob$bx <- rbind(blx=-b[(d+1):(2*d)], 
                   bux=b[1:d])
  
  for (j in 1:d) {
    
    prob$sense <- "min"
    prob$c <- A[j, ]
    r <- mosek(prob, list(verbose=0))
    stopifnot(identical(r$response$code, 0))
    min_dist = prob$c %*% r$sol$itr$xx
    #print(paste0("min_dist = ", min_dist))
    
    
    prob$sense <- "max"
    r <- mosek(prob, list(verbose=0))
    stopifnot(identical(r$response$code, 0))
    max_dist = prob$c %*% r$sol$itr$xx
    #print(paste0("mix_dist = ", max_dist))
    #print(" ")
    
    minFluxes = c(minFluxes, min_dist)
    maxFluxes = c(maxFluxes, max_dist)
    
    width = abs(max_dist - min_dist)
    
    if (width < 1e-07) {
      Aeq = rbind(Aeq, A[j,])
      beq = c(beq, min_dist)
    }
  }
  
  ret_list = list()
  ret_list$Aeq = Aeq
  ret_list$beq = beq
  ret_list$minFluxes = minFluxes
  ret_list$maxFluxes = maxFluxes
  
  return(ret_list)
}
