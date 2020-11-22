library(volesti)
library(Matrix)
library(Rmosek)

names = c('iAB_RBC_283', 'iAT_PLT_636', 'iAF1260', 'iBWG_1329', 'iEC1344_C', 'iJO1366', 'iML1515', 'recon1', 'recon2', 'recon3')

for (name in names) {
  

  load(paste0(name,'_polytope.RData'))

  d = dim(poly_list$A)[2]
  m = dim(poly_list$A)[1]
  N = 1000

  tim = system.time({ samples_list = multiphase_sampling(A = poly_list$A, b = poly_list$b, 
                                  max_ball = poly_list$max_ball, n = N,
                                  num_rounding_samples = 20*d, 
                                  max_num_samples = 100*d, rounding = TRUE) })
  
  N = dim(samples_list$samples)[2]
  if (poly_list$scale_applied) {
    samples = poly_list$T_scale %*% samples_list$samples + 
      kronecker(matrix(1, 1, N), matrix(poly_list$scale_shift, ncol = 1))
  } else {
    samples = samples_list$samples
  }

  steady_states = poly_list$N %*% samples + 
      kronecker(matrix(1, 1, N), matrix(poly_list$N_shift, ncol = 1))
  
  save(steady_states, file = paste0(name,'_steady_states.RData'))
  save(tim, file = paste0(name,'_sampling_time.RData'))
  #for (i in 1:N) {
  #  q = poly_list$A%*%samples_list$samples[, i] - poly_list$b
  #  if(length(which(q > 0))){
  #    print(i)
  #    #break
  #  }
  #}
  steady_states = matrix(0,0,0)
}


