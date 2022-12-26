mmc_sampling <- function(A, b, max_ball, n, num_rounding_samples, 
                                max_num_samples, rounding = TRUE) {
  
  d = dim(A)[2]
  m = dim(A)[1]
  
  parameters = list()
  parameters$A = A
  parameters$b = b
  parameters$T = diag(d)
  parameters$T_shift = rep(0, d)
  parameters$Neff = n
  parameters$Neff_sampled = 0
  parameters$num_rounding_steps = num_rounding_samples
  parameters$max_num_samples = max_num_samples
  parameters$round_it = 1
  parameters$request_rounding = rounding
  parameters$L = 6*sqrt(d)*max_ball$radius
  parameters$window = 100
  parameters$complete = FALSE
  parameters$rounding_completed = FALSE
  parameters$walk_length = 1
  
  samples = matrix(0, d, 0)
  T_curr = parameters$T
  T_shift_curr = parameters$T_shift
  neff_sampled = 0
  do_update = TRUE
  
  phase = 0
  while (!parameters$complete) {
    phase = phase + 1
    parameters = mmcs_phase(max_ball$center, max_ball$radius, parameters)
    #print(paste0("phase = ", phase,' performed'))
    N = dim(parameters$correlated_samples)[2]
    neff_sampled = neff_sampled + parameters$Neff_sampled
    
    parameters$Neff = parameters$Neff - parameters$Neff_sampled
    #print(paste0("Neff = ",  parameters$Neff ))
    
    if (parameters$request_rounding) {
      p = T_curr %*% parameters$correlated_samples[, 1:parameters$total_samples] + 
                      kronecker(matrix(1, 1, parameters$total_samples), matrix(T_shift_curr, ncol = 1))
      if (do_update) {
        T_curr = parameters$T
        T_shift_curr = parameters$T_shift
        max_ball = get_max_inner_ball_2(parameters$A, parameters$b)
      }
      if (parameters$rounding_completed) { 
        do_update = FALSE
      }
      parameters$L = 6*sqrt(d)*max_ball$radius
    } else {
      p = parameters$correlated_samples[, 1:parameters$total_samples]
    }
    samples = cbind(samples, p)
  }
  
  result_list = list()
  result_list$samples = samples
  result_list$Neff = neff_sampled
  result_list$A_rounded = parameters$A
  result_list$b_rounded = parameters$b
  result_list$T = T_curr
  result_list$T_shift = T_shift_curr
  
  return(result_list)
}