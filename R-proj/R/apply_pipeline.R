apply_pipeline <- function(path) {
  
  P = metabolic_net_2_polytope(path)
  print("compute min and max Fluxes")
  pre_proc_list = preprocess_with_mosek(P)
  
  print("Compute the null space to constraint")
  rr = full_dimensional_polytope_with_arma(pre_proc_list$Aeq, pre_proc_list$beq)
  
  A = P$A 
  b = P$b
  b = b - A %*% rr$N_shift
  A = A %*% rr$N
  
  m = dim(A)[1]
  d = dim(A)[2]
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
  
  HP = Hpolytope$new(A = A, b = b)
  ret_list = rounding_isotropic(HP)
  
  HP = Hpolytope$new(A = ret_list$A, b = ret_list$b)
  z=get_max_inner_ball(HP$A, HP$b)
  
  N = 3000
  print("Sample points from full diemnsional polytope")
  samples =  sample_points(HP, random_walk = list("walk" = "aBiW", "starting_point" = z$center,
                           "walk_length" = 1, "L" = 2*d*z$radius), n = N)
  
  samples = ret_list$T %*% samples + 
                  kronecker(matrix(1, 1, N), matrix(ret_list$T_shift, ncol = 1))
  
  steady_states = rr$N %*% samples + 
                  kronecker(matrix(1, 1, N), matrix(rr$N_shift, ncol = 1))
  
  result_list = list()
  result_list$HP_rounded = HP
  result_list$steady_states = steady_states
  result_list$N = rr$N
  result_list$N_shift = rr$N_shift
  result_list$T = ret_list$T
  result_list$T_shift = ret_list$T_shift
  
  return(result_list)
}
