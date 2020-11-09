preprocess_facet_widths <- function(P) {
  
  d = P$dimension
  Aeq = P$Aeq
  beq = P$beq
  
  A = P$A
  b = P$b
  
  lps.model <- lpSolveAPI::make.lp(0, P_low$dimension)
  
  for (i in 1:dim(Aeq)[1]) {
    lpSolveAPI::add.constraint(lps.model, Aeq[i,], "=", beq[i])
  }
  lpSolveAPI::set.bounds(lps.model, lower = -b[(d+1):(2*d)], upper = b[1:d], columns = 1:d)
  
  Aeq_new = Aeq
  beq_new = beq
  
  A_new = matrix(0,0,d)
  b_new = c()
  
  for (j in 1:dim(A)[1]) {
    
    lpSolveAPI::set.objfn(lps.model, A[j,])
    solve(lps.model)
    max_dist = lpSolveAPI::get.objective(lps.model)
    lpSolveAPI::set.objfn(lps.model, -A[j,])
    solve(lps.model)
    min_dist = lpSolveAPI::get.objective(lps.model)
    width = abs(max_dist+min_dist)/norm(as.matrix(A[i,]));
    
    if (width<1e-07) {
      Aeq_new = rbind(Aeq_new, A[j,])
      beq_new = c(beq_new, max_dist)
    } else {
      A_new = rbind(A_new, A[j,])
      b_new =c(b_new, b[j])
    }
  }
  
  PP = Hpolytope$new(A=A_new, b=b_new, Aeq = Aeq_new, beq = beq_new)
  return (PP)
}


