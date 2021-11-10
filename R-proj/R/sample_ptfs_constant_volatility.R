#' @export
sample_ptfs_constant_volatility <- function(sigma, c, M, parameters) {
  
  nu = parameters$nu
  lb = parameters$lb
  ub = parameters$ub
  Nu = parameters$Nu
  W = parameters$W
  
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

  res = get_tree_of_components(V, A, b, center_2, lb, ub, nu, Nu, W)
  
  tree = res$tree
  single_node = res$single_node
  
  if( single_node) {
    print('height of tree')
    #print(count_height_of_tree(tree))
  }
  samples = get_samples_from_tree(tree, single_node, A, b, center_2, V, M, parameters$W_to_sample, parameters$psrf_target, parameters$error, parameters$win_len)
                                 #(tree, single_node, A, b, x0, V, N, W, psrf_target)
  
  NN = dim(samples)[2]
  
  samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
  #print(dim(samples))
  #print(dim(N))
  #print(length(center))
  #print(n)
  #print(dim(kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))))
  samples = N %*% (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
  #print(dim(samples))
  #samples = samples 
  
  return(samples)
}
  