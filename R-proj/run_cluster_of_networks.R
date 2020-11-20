names = c()

for (name in names) {
  

  load(paste0(name,'_polytope.RData'))

  d = dim(A)[2]
  m = dim(A)[1]
  N = 1000

  tim = system.time({ samples = multiphase_sampling(A = A, b = b, 
                                  max_ball = max_ball, n = N,
                                  num_rounding_samples = 10*d, 
                                  max_num_samples = 10*d, rounding = TRUE) })

  steady_states = N %*% samples + 
      kronecker(matrix(1, 1, N), matrix(N_shift, ncol = 1))
  
  save(steady_states, file = paste0(name,'_steady_states.RData'))
  save(tim, file = paste0(name,'_sampling_time.RData'))
  steady_states = matrix(0,0,0)
}