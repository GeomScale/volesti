apply_pipeline <- function(path) {
  
  P = metabolic_net_2_polytope(path)
  pre_proc_list = preprocess_with_mosek(P)
  rr = full_dimensional_polytope_with_arma(pre_proc_list$Aeq, pre_proc_list$beq)
  
  A = P$A 
  b = P$b
  b = b - A%*%rr$N_shift
  A = A %*% rr$N
  
  m = dim(A)[1]
  d = dim(A)
  rows_to_del = c()
  for (i in 1:m) {
    if (sqrt(sum(A[i,]^2)) < 1e-06) {
      rows_to_del = c(rows_to_del, i)
    }
  }
  if (length(rows_to_del) > 0) {
    A = A[-rows_to_del, ]
    b = b[-rows_to_del]
  }
  
  z=get_max_inner_ball(A, b)
  return(z)
}