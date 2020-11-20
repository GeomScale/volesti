library(volesti)
library(Matrix)
library(Rmosek)

path_mets = "/home/tolis/data/metabolic_mat/mat_socg/e_coli"
file.names <- dir(path_mets, pattern =".mat")
source('~/volume_approximation/R-proj/R/null_space_and_shift.R', echo=TRUE)
time_vec = c()

for (i in 1:length(file.names)) {
  
  path = paste0(path_mets, "/", file.names[i])
  print(path)
  name = file.names[i]
  name = gsub('.{4}$','',name)
  print(name)
  
  tim = system.time({
  P = metabolic_net_2_polytope(path)
  print("compute min and max Fluxes")
  pre_proc_list = fast_preprocess_with_mosek(P)
  save(pre_proc_list, file = paste0(name,'_pre_process.RData'))
  
  print("Compute the null space to constraint")
  #rr = full_dimensional_polytope_with_arma(pre_proc_list$Aeq, pre_proc_list$beq)
  rr = null_space_and_shift(pre_proc_list$row_ind, pre_proc_list$col_ind, pre_proc_list$values, pre_proc_list$Aeq, pre_proc_list$beq)
  print("null space computed")
  save(rr, file = paste0(name,'_null_space.RData'))
  
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
  print(paste0(length(rows_to_del), " facets removed"))
  
  print("Compute scaling for numerical stability")
  sc = central_scaling(A, b) 
  })
  
  A = sc$A
  b = sc$b
  T_scale = sc$T_scale
  scale_shift = sc$scale_shift
  save(A, b, sc, file = paste0(name,'_A_b_ball.RData'))
  save(tim, file = paste0(name,'_total_preprocess_time.RData'))
  
  max_ball = list()
  max_ball$center = sc$center
  max_ball$radius = sc$radius
  
  poly_list = list()
  poly_list$A = A
  poly_list$b = b
  poly_list$N = rr$N
  poly_list$N_shift = rr$N_shift
  poly_list$max_ball = max_ball
  poly_list$scale_applied = sc$scale_applied 
  poly_list$T_scale = sc$T_scale
  poly_list$scale_shift = sc$scale_shift
  
  save(poly_list, file = paste0(name,'_polytope.RData'))
  
}
