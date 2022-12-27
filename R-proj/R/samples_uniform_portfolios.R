samples_uniform_portfolios <- function(A, b, Aeq, beq, ess = 1000) {
  
  print("Preprocessing...")
  pre_proc_list = preprocess_with_quadprog(A, b, Aeq, beq)
  
  print("Computing the full dimensional polytope...")
  rr = null_space_and_shift(pre_proc_list$Aeq, pre_proc_list$beq)
  
  b = b - A %*% rr$N_shift
  A = A %*% rr$N
  
  m = dim(A)[1]
  d = dim(A)[2]

  print('Full dimensional polytope computed!')
  
  print("Computing the Chebychev ball of the full dimensional polytope...")
  max_ball = get_max_inner_ball(A, b)
  print("Chebychev ball computed!")
  
  d = dim(A)[2]
  m = dim(A)[1]
  
  print("Sampling with MMCS method...")
  tim = system.time({ samples_list = mmc_sampling(A = A, b = b, 
                                                  max_ball = max_ball, n = ess,
                                                  num_rounding_samples = 20*d, 
                                                  max_num_samples = 100*d, 
                                                  rounding = TRUE) })
  print('Sampling completed!')
  
  N = dim(samples_list$samples)[2]
  samples = samples_list$samples
  
  steady_states = rr$N %*% samples + 
    kronecker(matrix(1, 1, N), matrix(rr$N_shift, ncol = 1))
  
  HP = Hpolytope$new(A = samples_list$A_rounded, b =  samples_list$b_rounded)
  
  result_list = list()
  result_list$HP_rounded = HP
  result_list$samples = samples
  result_list$random_portfolios = steady_states
  result_list$N = rr$N
  result_list$N_shift = rr$N_shift
  result_list$T = samples_list$T
  result_list$T_shift = samples_list$T_shift
  result_list$minWeights = pre_proc_list$minWeights
  result_list$maxWeights = pre_proc_list$maxWeights
  result_list$run_time = tim[3]
  
  return(result_list)
}
