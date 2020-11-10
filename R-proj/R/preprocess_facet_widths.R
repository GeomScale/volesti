preprocess_facet_widths <- function(P) {
  
  d = P$dimension
  Aeq = P$Aeq
  beq = P$beq
  
  A = P$A
  b = P$b
  
  lps.model <- lpSolveAPI::make.lp(0, P$dimension)
  
  for (i in 1:dim(Aeq)[1]) {
    lpSolveAPI::add.constraint(lps.model, Aeq[i,], "=", beq[i])
  }
  lpSolveAPI::set.bounds(lps.model, lower = -b[(d+1):(2*d)], upper = b[1:d], columns = 1:d)
  
  Aeq_new = Aeq
  beq_new = beq
  
  A_new = matrix(0,0,d)
  b_new = c()
  
  minFluxes = c()
  maxFluxes = c()
  
  for (j in 1:d) {
    
    lpSolveAPI::set.objfn(lps.model, A[j,])
    solve(lps.model)
    min_dist = lpSolveAPI::get.objective(lps.model)
    lpSolveAPI::set.objfn(lps.model, -A[j,])
    solve(lps.model)
    max_dist = -lpSolveAPI::get.objective(lps.model)
    width = abs(max_dist-min_dist)
    
    minFluxes = c(minFluxes, min_dist)
    maxFluxes = c(maxFluxes, max_dist)
    
    if (width<1e-07) {
      Aeq_new = rbind(Aeq_new, A[j,])
      beq_new = c(beq_new, min_dist)
    }
  }
  
  ret_list = list()
  ret_list$Aeq = Aeq_new
  ret_list$beq = beq_new
  ret_list$minFluxes = minFluxes
  ret_list$maxFluxes = maxFluxes
  
  return(ret_list)
}


