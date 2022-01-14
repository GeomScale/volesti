#' @export
sample_ptfs_constant_volatility <- function(sigma, c, M, ignore_smallest_components = TRUE) {
  
  n <- dim(sigma)[2]
  
  parameters = get_parameters(n-1)
  
  nu = parameters$nu
  lb = parameters$lb
  ub = parameters$ub
  Nu = parameters$Nu
  W = parameters$W
  psrf_target = parameters$psrf_target
  error = parameters$error
  WW = parameters$WW
  
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
  
  components_res=find_components_new_vertices(V, center_2)
  
  S = components_res[[1]]
  S_Vindices = components_res[[2]]
  
  if(length(components_res) > 2){
    V_ind_out = components_res[[3]]
  } else {
    V_ind_out = c()
  }
  
  cmax = get_c_upper_bound(A, b, center_2)
  cmin = 1
  
  print(S_Vindices)
  #print(W)
  
  res = list()
  res$overall_samples = matrix(list(), 0, 1)
  res$num_verts = matrix(list(), 0, 1)
  
  if (length(S) == 1) {
    if (length(S_Vindices[[1]]) == dim(V)[2]) {
      interior_res = compute_interior_point_in_node_eff(S_Vindices[[1]], V_ind_out, V, center_2, A, b)
      if (interior_res$found) {
        x = interior_res$x
      } else {
        x = compute_interior_point_single_component(A, b, center_2)
      }
      #samples = sample_component_psrf(A, b, x, 2000, W, center_2, psrf_target)
      samples = sample_component_psrf_interface(A, b, x, M, W, center_2, psrf_target, V, S_Vindices[[1]])
      
      if (dim(samples)[2] > M) {
        indx <- sample(1:dim(samples)[2], M)
        samples = samples[, indx]
      }
      NN = dim(samples)[2]
      samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
      samples = N %*% (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
      
      res$overall_samples[[length(res$overall_samples) + 1]] = samples
      res$num_verts[[length(res$num_verts) + 1]] = length(S_Vindices[[1]])
      
      return(res)
      
    } else {
      
      print("sample from one component")
      Xs = get_points_on_component_imp(A, b, center_2, V, S_Vindices[[1]])
      nnn = dim(Xs)[2]
      print(paste0("number of starting points = ",as.character(nnn)))
      d = length(center_2)
      samples = matrix(,d,0)

      for (ii in 1:nnn) {
        X = sample_component_psrf_interface(A, b, Xs[,ii], M, W, center_2, psrf_target, V, S_Vindices[[1]])
        if (dim(X)[2] > M) {
          indxes <- sample(1:dim(X)[2], M)
          X = X[, indxes]
        }
        samples = cbind(samples, X)
      }
      
      MM = 50000
      if (dim(samples)[2] > MM) {
        indxes <- sample(1:dim(samples)[2], MM)
        samples = samples[, indxes]
      }
      
      NN = dim(samples)[2]
      #print(paste0("num_samples = ",as.character(NN),", num_vertices = ", as.character(length(S_Vindices[[indx]])),", Mxnv = ",M*length(S_Vindices[[indx]])))
      samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
      samples = N %*% (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
      
      res$overall_samples[[length(res$overall_samples) + 1]] = samples
      res$num_verts[[length(res$num_verts) + 1]] = length(S_Vindices[[1]])
      
      return(res)
    }
  }
  
  nn = length(S)
  max_num = 0
  index_max = 0
  
  for (jj in 1:nn) {
    Sind = S_Vindices[[jj]]
    if (length(Sind) > max_num){
      index_max = jj
      max_num = length(Sind)
    }
  }
  
  print(paste0("max_num_verts = ",as.character(max_num), ", index_num = ",as.character(index_max)))
  
  Sind = S_Vindices[[index_max]]
  xc = get_center_max_ball(A, b, center_2, V, Sind)
  
  for (i in 1:nn) {
    
    if(ignore_smallest_components) {
      #print(paste0("number of vertices = ",as.character(length(S_Vindices[[i]]))))
      if (i != index_max) {
        #print("component ignored")
        next
      }
    }
    
    print(paste0("index of component = ",as.character(i)))
    Xs = get_points_on_component_imp(A, b, center_2, V, S_Vindices[[i]], xc)
    nnn = dim(Xs)[2]
    print(paste0("number of starting points = ",as.character(nnn)))
    
    d = length(center_2)
    samples = matrix(,d,0)

    for (ii in 1:nnn) {
      X = sample_component_psrf_interface(A, b, Xs[,ii], M, W, center_2, psrf_target, V, S_Vindices[[i]])
      if (dim(X)[2] > M) {
        indxes <- sample(1:dim(X)[2], M)
        X = X[, indxes]
      }
      samples = cbind(samples, X)
    }
    
    MM = 5000
    if (nnn > 1) {
      MM = 50000
    }
    
    if (dim(samples)[2] > MM) {
      indxes <- sample(1:dim(samples)[2], MM)
      samples = samples[, indxes]
    }
    
    NN = dim(samples)[2]
    #print(paste0("num_samples = ",as.character(NN),", num_vertices = ", as.character(length(S_Vindices[[indx]])),", Mxnv = ",M*length(S_Vindices[[indx]])))
    samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
    samples = N %*% (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
    
    res$overall_samples[[length(res$overall_samples) + 1]] = samples
    res$num_verts[[length(res$num_verts) + 1]] = length(S_Vindices[[i]])
  }
  
  return(res)
}

