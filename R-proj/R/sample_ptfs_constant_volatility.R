#' @export
sample_ptfs_constant_volatility <- function(sigma, c, parameters) {
  
  nu = parameters$nu
  lb = parameters$lb
  ub = parameters$ub
  Nu = parameters$Nu
  W = parameters$W
  
  n <- (dim(sigma)[2])
  
  A <- -diag(n)
  b <- rep(0, 1)
  
  Aeq <- rep(1, n)
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
  
  R = c + center%*%sigma_proj%*%center - (rep(1,n)/n)%*%sigma%*%(rep(1,n)/n)
  sigma_proj = sigma_proj / R
  
  b = b - A %*% center # shift ellipsoid's center to origin
  x0 = x0 - center
  V = V - kronecker(matrix(1, 1, n), matrix(center, ncol = 1))
  center_2 = rep(0, n-1)
  
  Tinv = t(chol(MASS::ginv(sigma_proj)))
  T = MASS::ginv(Tinv)

  A = A %*% Tinv
  x0_2 = T %*% x0
  V = T %*% V
    
  b = b - A %*% x0_2
  center_2 = center_2 - x0_2
  V = V - kronecker(matrix(1, 1, n), matrix(x0_2, ncol = 1))

  res = get_tree_of_components(V, A, b, center_2, lb, ub, nu, Nu, W)
  
  tree = res$tree
  single_node = res$single_node
  
  samples = get_samples_from_tree(tree, single_node, A, b, center_2, Nu, W)
  
  NN = dim(samples)[2]
  samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
  samples = N * (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
  
  return(samples)
}
  