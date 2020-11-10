apply_pipeline <- function(path, remove_biomass = FALSE) {
  
  P = metabolic_net_2_polytope(path, remove_biomass)
  print("compute min and max Fluxes")
  pre_proc_list = fast_preprocess_with_mosek(P)
  
  print("Compute the null space to constraint")
  rr = full_dimensional_polytope_with_arma(pre_proc_list$Aeq, pre_proc_list$beq)
  
  print("Get full dimensional polytope")
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
  A = sc$A
  b = sc$b
  T_scale = sc$T_scale
  scale_shift = sc$scale_shift
  
  #z=get_max_inner_ball(A, b)
  #reduce_factor = 1
  #done = FALSE
  #max_iter = ceiling(log10(1/(z$radius[1,1]))) - 1
  #iter = 0
  #T_scale = diag(d)
  #scale_shift = rep(0, d)
  #while(!done & max_iter > 0) {
  #  iter = iter + 1
  #  print(paste0("iter = ",iter,", max_iter = ", max_iter))
  #  scale_shift = z$center
  #  brep = b - A %*% scale_shift
  #  scale_factor = 1/(z$radius[1,1]) 
  #  scale_factor = scale_factor / reduce_factor
  #  brep2 = brep * (scale_factor)
  #  a = try(get_max_inner_ball(A, brep2))
  #  if (class(a)=="list") {
  #    print(a$radius)
  #    print(scale_factor)
  #    if (a$radius > z$radius) {
  #      done = TRUE
  #      b = brep2
  #      T_scale = diag(d) * scale_factor
  #    } else {
  #      reduce_factor = reduce_factor * 10
  #    }
  #  } else {
  #    reduce_factor = reduce_factor * 10
  #  }
  #  if (iter == max_iter) {
  #    T_scale = diag(d)
  #    scale_shift = rep(0, d)
  #    print("we exceed maximmum number of iterations for scaling")
  #    break
  #  }
  #}
  ###----scaling ends here----##
  
  print("Rounding the polytope")
  HP = Hpolytope$new(A = A, b = b)
  ret_list = rounding_isotropic(HP)
  
  HP = Hpolytope$new(A = ret_list$A, b = ret_list$b)
  z=get_max_inner_ball(HP$A, HP$b)
  
  N = 3000
  print("Sample points from full diemnsional polytope")
  samples =  sample_points(HP, random_walk = list("walk" = "aBiW", "starting_point" = z$center,
                           "walk_length" = 1, "L" = 4*sqrt(d)*z$radius), n = N)
  
  samples = ret_list$T %*% samples + 
                  kronecker(matrix(1, 1, N), matrix(ret_list$T_shift, ncol = 1))
  
  samples = T_scale %*% samples + 
                  kronecker(matrix(1, 1, N), matrix(scale_shift, ncol = 1))
  
  steady_states = rr$N %*% samples + 
                  kronecker(matrix(1, 1, N), matrix(rr$N_shift, ncol = 1))
  
  result_list = list()
  result_list$HP_rounded = HP
  result_list$steady_states = steady_states
  result_list$N = rr$N
  result_list$N_shift = rr$N_shift
  result_list$T = ret_list$T
  result_list$T_shift = ret_list$T_shift
  result_list$minFluxes = pre_proc_list$minFluxes
  result_list$maxFluxes = pre_proc_list$maxFluxes
  
  return(result_list)
}
