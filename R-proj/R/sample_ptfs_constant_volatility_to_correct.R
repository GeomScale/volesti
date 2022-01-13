#' @export
sample_ptfs_constant_volatility_to_correct <- function(sigma, c, M) {
  
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
  print(W)
  
  res = list()
  
  if (length(S) == 1 && length(S_Vindices[[1]]) != dim(V)[2]){
    
    Xs = get_points_on_components_imp(A, b, center_2, V, S, S_Vindices)
    
    nn = length(S)
    a_vals = matrix(list(), 0, 1)
    ratios =  matrix(list(), 0, 1)
    mus = matrix(list(), 0, 1)
    sigmas = matrix(list(), 0, 1)
    XX = matrix(list(), 0, 1)
    vols = c()
    num_verts=c()
    
    for (i in 1:nn) {
      num_verts = c(num_verts, length(S_Vindices[[i]]))
    }
    num_verts2 = sort(num_verts)
    tails_v = tail(num_verts2, 2)
    #print(paste0('num_verts = ',as.character(num_verts)))
    #print(paste0('num_verts2 = ',as.character(num_verts2)))
    print(paste0('tails_v = ',as.character(tails_v)))
    
    indices = which(num_verts == max(num_verts))
    #print(indices)
    indx = which(num_verts == max(num_verts))
    for (i in 1:nn) {
      mu = Xs[,i]
      mus[[length(mus) + 1]] = mu 
    }
    
    print("sample from one component")
    print(paste0("indx = ",as.character(indx)))
    Xs = get_points_on_component_imp(A, b, center_2, V, S_Vindices[[indx]])
    nnn = dim(Xs)[2]
    print(paste0("number of starting points = ",as.character(nnn)))
    d = length(center_2)
    samples = matrix(,d,0)
    for (ii in 1:nnn) {
      X = sample_component_psrf_interface(A, b, Xs[,ii], M, W, center_2, psrf_target, V, S_Vindices[[indx]])
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
    print(paste0("num_samples = ",as.character(NN),", num_vertices = ", as.character(length(S_Vindices[[indx]])),", Mxnv = ",M*length(S_Vindices[[indx]])))
    samples = Tinv %*% (samples - kronecker(matrix(1, 1, NN), matrix(center_2, ncol = 1)))
    samples = N %*% (samples + kronecker(matrix(1, 1, NN), matrix(center, ncol = 1))) + kronecker(matrix(1, 1, NN), matrix(rep(1,n)/n, ncol = 1))
    
    res$computations = TRUE
    res$samples = samples
    
    return(res)
  }
  
  
  
  res$computations = FALSE
  res$samples = c()
  return(res)
  
  
}


