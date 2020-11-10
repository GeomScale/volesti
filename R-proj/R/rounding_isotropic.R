rounding_isotropic <- function(P) {
  
  A = P$A
  b = P$b
  
  d = dim(A)[2]
  m = dim(A)[1]
  
  parameters = list()
  parameters$A = A
  parameters$b = b
  parameters$T = diag(d)
  parameters$T_shift = rep(0, d)
  parameters$num_rounding_steps = 10 * d
  parameters$round_it = 1
  parameters$max_s = 1e+50
  parameters$prev_max_s = 1e+50
  parameters$fail = FALSE
  parameters$converged = FALSE
  parameters$last_round_under_p = FALSE
  walk_length = 2
  
  while (!parameters$converged) {
    ball = get_max_inner_ball(parameters$A, parameters$b)
    center = ball$center
    r = ball$radius
    parameters = rounding_svd_step(center, r, walk_length, parameters)
  }
  
  result_list = list()
  result_list$A = parameters$A
  result_list$b = parameters$b
  result_list$T = parameters$T
  result_list$T_shift = parameters$T_shift
  
  return(result_list)
}