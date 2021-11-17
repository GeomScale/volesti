#' @export
sample_ptfs_constant_volatility <- function(sigma, c, M, parameters) {
  
  nu = parameters$nu
  lb = parameters$lb
  ub = parameters$ub
  Nu = parameters$Nu
  W = parameters$W
  psrf_target = parameters$psrf_target
  error = parameters$error
  WW = parameters$WW
  
  n <- dim(sigma)[2]
  
  A <- -diag(n)
  b <- rep(0, 1)
  
  Aeq <- matrix(rep(1, n), nrow = 1, ncol = n)
  beq <- c(1)
  
  N = pracma::nullspace(Aeq)
  x0 <- rep(1, n)/n
  
  b = b - A %*% x0
  A = A %*% N
  V = diag(n)
  V = V - kronecker(matrix(1, 1, n), matrix(x0, ncol = 1))
  V = t(N) %*% V
  x0 = rep(0, n-1)
  
  sigma_proj = t(N) %*% sigma %*% N
  center = -MASS::ginv(sigma_proj) %*% (t(N) %*% sigma) %*% (rep(1,n)/n)
  
  R = c + t(center)%*%sigma_proj%*%center - (rep(1,n)/n)%*%sigma%*%(rep(1,n)/n)
  R = R[1]
  sigma_proj = sigma_proj / R
  
  b = b - A %*% center # shift ellipsoid's center to origin
  x0 = x0 - center
  V = V - kronecker(matrix(1, 1, n), matrix(center, ncol = 1))
  center_2 = rep(0, n-1)
  
  Tinv = t(chol(Matrix::nearPD(MASS::ginv(sigma_proj))$mat))
  T = MASS::ginv(Tinv)
  
  A = A %*% Tinv
  x0_2 = T %*% x0
  V = T %*% V
  
  b = b - A %*% x0_2
  center_2 = center_2 - x0_2
  V = V - kronecker(matrix(1, 1, n), matrix(x0_2, ncol = 1))
  
  find_components_new_vertices(V, center_2)
  
  S = components_res[[1]]
  S_Vindices = components_res[[2]]
  
  if(length(components_res) > 2){
    V_ind_out = components_res[[3]]
  } else {
    V_ind_out = c()
  }
  
  if (length(S) == 1){
    interior_res = compute_interior_point_in_node_eff(S_Vindices[[1]], V_ind_out, V, center_2, A, b)
    if (interior_res$found) {
      x = interior_res$x
    } else {
      x = get_fast_interior_point(A, b, center_2, V, S_Vindices[[1]], cmin, cmax)
    }
    samples = sample_component_psrf(A, b, x, 2*Nu, W, center_2, psrf_target)
    
    if (dim(samples)[2] > Nu) {
      indx <- sample(1:dim(samples)[2], Nu)
      samples = samples[,indx]
    }
    NN = dim(samples)[2]
    samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
    samples = N %*% (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
    return(samples)
  }
  
  Xs = get_points_on_components(A, b, center_2, V, S, S_Vindices, V_ind_out)
  
  nn = length(S)
  a_vals_max = c()
  ratios = c()
  mus = c()
  XX = matrix(list(), 0, 1)
  
  for (i in 1:nn) {
    X = sample_component_psrf(A, b, Xs[,i], 2*Nu, W, center_2, psrf_target)
    XX[[length(XX)+1]] = X
    
    if (50*n < dim(X)[2]){
      NN = 50*n
    } else {
      NN = dim(X)[2]
    }
    
    res = get_max_cap(X, A, b, center_2, V, S_Vindices[[i]], NN)
    mu = res$center
    sigma = cov(t(X))
    
    res = compute_component_volume_fischer(A, b, mu, sigma, center_2, X, WW, error)
    a_vals_max = c(a_vals_max, res$last_variance)
    ratios = c(ratios, res$ratio)
  }
  
  relative_vols = get_rel_volumes(a_vals_max, ratios, A, b, mu, sigma)
  samples = get_samples(XX, relative_vols, Nu)
  
  NN = dim(samples)[2]
  samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
  samples = N %*% (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
  
  return(samples)
}
