generate_steady_states <- function(path, n = 1000, Recon2D_v04 = FALSE, Recon3D_301 = FALSE) {
  
  P = metabolic_net_2_polytope(path, Recon2D_v04, Recon3D_301)
  print("computing min and max Fluxes")
  pre_proc_list = fast_preprocess_with_mosek(P)
  
  print("Computing the null space of the Stoichiometric matrix")
  rr = null_space_and_shift(pre_proc_list$row_ind, pre_proc_list$col_ind, pre_proc_list$values, pre_proc_list$Aeq, pre_proc_list$beq)
  
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
  print('Full dimensional polytope computed!')
  
  print("Computing the Chebychev ball of the full dimensional polytipe")
  max_ball = get_max_inner_ball(A, b)
  
  d = dim(A)[2]
  m = dim(A)[1]
  
  print("Call MMCS method")
  tim = system.time({ samples_list = mmc_sampling(A = A, b = b, 
                                                         max_ball = max_ball, n = n,
                                                         num_rounding_samples = 20*d, 
                                                         max_num_samples = 100*d, 
                                                         rounding = TRUE) })
  print('Steady states computed!')
  
  N = dim(samples_list$samples)[2]
  samples = samples_list$samples
  
  steady_states = rr$N %*% samples + 
    kronecker(matrix(1, 1, N), matrix(rr$N_shift, ncol = 1))
  
  HP = Hpolytope$new(A = samples_list$A_rounded, b =  samples_list$b_rounded)
  
  result_list = list()
  result_list$HP_rounded = HP
  result_list$samples = samples
  result_list$steady_states = steady_states
  result_list$N = rr$N
  result_list$N_shift = rr$N_shift
  result_list$T = samples_list$T
  result_list$T_shift = samples_list$T_shift
  result_list$minFluxes = pre_proc_list$minFluxes
  result_list$maxFluxes = pre_proc_list$maxFluxes
  
  return(result_list)
}
