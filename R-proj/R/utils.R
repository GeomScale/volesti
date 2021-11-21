
#' @export
min_variance <- function(Dmat) {
  
  Amat <- matrix(0, nrow = ncol(Dmat) * 2 + 1, ncol = ncol(Dmat),
                 dimnames = list(NULL, colnames(Dmat)) )
  Amat[1, ] <- 1
  Amat[2:(ncol(Dmat)+1), ] <- diag(ncol(Dmat))
  Amat[(ncol(Dmat)+2):nrow(Amat), ] <- -diag(ncol(Dmat))
  Amat <- t(Amat)
  bvec <- c(1, rep(0, ncol(Dmat)), rep(-1, ncol(Dmat)) )
  dvec <- rep(0, ncol(Dmat))
  
  opt <- quadprog::solve.QP( Dmat = Dmat,
                             dvec = dvec,
                             Amat = Amat,
                             bvec = bvec,
                             meq  = 1 )
  
  x = opt$solution
  
  min_var = t(x) %*% sigma %*% x
  
  return(min_var[1])
}


#' @export
max_variance <- function(sigma) {
  
  diag_sigma = diag(sigma)
  pos = which(diag_sigma == max(diag_sigma))
  x = rep(0, ncol(sigma))
  x[pos] = 1
  max_var = t(x) %*% sigma %*% x
  
  return(max_var[1])
}



get_sequence_of_volatilities_CB <- function( sigma, m, ignore_cov = TRUE )
{
  if ( isTRUE(ignore_cov) ) {
    # Override covariance coefficients with zero to ensure that
    # the output (m portfolio variances) are monotonely increasing

    sigma[upper.tri(sigma)] <- 0
    sigma[lower.tri(sigma)] <- 0
  }
  
  # Sort by increasing variances
  sds <- diag( sigma )
  ordering <- order(sds)
  
  # Estimate the variance of an equally-weighted portfolio within
  # each vola-quantile
  
  portfolio_variances <- numeric(m)
  
  for ( i in 1:m ) {
    k <- ceiling( length(ordering) / m )
    
    if ( i < m ) {
      idx <- ordering[(i*k-k+1):(i*k)]
    } else {
      idx <- ordering[(i*k-k+1):length(ordering)]
    }
    
    wghts <- rep(1/length(idx), length(idx))
    portfolio_variances[i] <- t(wghts) %*% sigma[idx, idx] %*% wghts
  }
  
  return( portfolio_variances )
}


#' @export
get_sequence_of_volatilities_2 <- function(sigma, m) {
  
  min_var = min_variance(sigma)
  max_var = max_variance(sigma)
  
  step = (max_var - min_var)/(m+2)
  Cs = seq(from=min_var, to = max_var, by=step)
  Cs = Cs[2:(length(Cs)-1)]
  
  return(Cs)
}


#' @export
get_sequence_of_volatilities <- function(sigma, m) {
  
  n = dim(sigma)[2]
  
  NN1 = 55000
  NN2 = 50000
  X = Sampling_simplex(n, NN1)
  
  Sx = sigma%*%X
  xSx = rep(0,NN1)
  
  for (j in 1:NN1) {
    xSx[j] = (t(X[,j]) %*% Sx[,j])[1]
  }
  I2 = order(xSx)
  Cs = rep(0,m)
  
  for (j in 1:m) {
    range = j*floor(NN2/m)
    Cs[j] = xSx[I2[range]]
  }
  return(Cs)
}


#' @export
check_correctness <- function(samples, sigma, c) {
  correctness = TRUE
  q = colSums(samples)
  if (length(q[which(q>1+1e-07 || q<1-1e-07)]) > 0) {
    print(length(q[which(q>1+1e-07 || q<1-1e-07)]))
    correctness = FALSE
  }
  
  N = dim(samples)[2]
  
  for (i in 1:N) {
    vol = (samples[,100] %*%sigma%*%samples[,100])[1]
    if (vol < (c - 1e-07) || vol > (c + 1e-07)){
      print(i)
      correctness = FALSE
    }
  }
  return(correctness)
}


#' @export
get_parameters <- function(d) {
  
  parameters = list()
  parameters$nu = 10
  parameters$lb = 0.1
  parameters$ub = 0.15
  parameters$Nu = 1000 + floor(d^2/2)
  parameters$W = 1
  parameters$WW = 7
  #parameters$W_to_sample = 10 + floor(d/10)
  parameters$W_to_sample = 1
  parameters$win_len = 4*d^2 + 500
  parameters$error = 0.1
  parameters$psrf_target = 1.2
  
  return(parameters)
  
}


#' @export
get_c_upper_bound <- function(A, b, x0) {

  row_norms = sqrt(rowSums(A^2))
  A = diag(1 / row_norms) %*% A
  b = diag(1 / row_norms) %*% b

  cmax = (A %*% x0 + 1) / b
  cmax = max(cmax) * 1.1

  return(cmax)
}

#' @export
compute_interior_point_in_node_eff <- function(S_Vindices, V_ind_out, V, x0, A, b) {

  found = TRUE
  res = list()

  q = A %*% x0 - b
  if (sum(q>0) == 0){
    x = V[, S_Vindices[1]]

    ball_line_res = ball_line_intersection(x, x0-x, x0, 1)
    t = c(ball_line_res$tmin, ball_line_res$tmax)
    t = min(t[t>0])
    p = x + t * (x0-x)
    res$found = found
    res$x = p
    return(res)
  }
  n = length(x0)

  if (length(V_ind_out) == 0) {
  
    p=c()
    found = FALSE
    res$found = found
    res$x = p
    return(res)
  }

  x1 = V[, S_Vindices[1]]
  x2 = V[, V_ind_out[1]]

  centroid = rep(0,n)#  zeros(n,1);
  if (sqrt(sum((centroid - x0)^2)) < 1) {
  
    v = centroid - x1

    #[tmin, tmax, ~] = ball_line_intersection(x1, v, x0, 1);
    ball_line_res = ball_line_intersection(x1, v, x0, 1)
    t = c(ball_line_res$tmin, ball_line_res$tmax)
    t = min(t[t>0])
    #t = [tmin; tmax];
    #t = min(t(t>0));
    p = x1 + t * v
  } else {
    v = centroid - x1

    ball_line_res = ball_line_intersection(x1, v, x0, 1)
    if (!ball_line_res$intersect || ((ball_line_res$tmin>1 || ball_line_res$tmin<0) & (ball_line_res$tmax>1 || ball_line_res$tmax<0))) {
  
      v = centroid - x2

      ball_line_res = ball_line_intersection(x2, v, x0, 1)
      t = c(ball_line_res$tmin, ball_line_res$tmax)
      t = min(t[t>0])
      p = x2 + t * v
    } else {
      v = centroid - x1

      ball_line_res = ball_line_intersection(x1, v, x0, 1)
      t = c(ball_line_res$tmin, ball_line_res$tmax)
      t = min(t[t>0])
      p = x1 + t * v
    }
  }
  
  res$found = found
  res$x = p
  return(res)
}


#' @export
GenerateNode <- function(c, x, V, indices, ratio) {
  
  node = list()
  node$c = c
  node$x0 = x
  node$V = V
  node$V_indices = indices
  node$ratio = ratio
  node$isLeaf = FALSE
  node$number_of_leaves = 0
  node$ratio_of_leaves = c()
  node$S_ind_to_esti = c()
  node$S_Vindices_leaves = matrix(list(), 0, 1)
  node$children = matrix(list(), 0, 1)
  node$error_to_leaf = NaN
  node$error_to_estimate = NaN
  node$max_depth = 0
  node$depth = 0
  
  return(node)
  
}

#' @export
boundary_randsphere <- function(d, N) {
  
  X = matrix( rnorm(d*N,mean=0,sd=1), d, N)
  a = sqrt(colSums(X^2))
  X = X / matrix(rep(a, d), ncol=N, byrow=T)
  
  return(X)
}


#' @export
remove_leaf <- function(S, S_Vindices, S_ind_to_esti) {
  
  n = length(S)
  S_new = matrix(list(), 0, 1)
  S_Vindices_new = matrix(list(), 0, 1)
  
  counter = 1
  for (i in 1:n) {
    if (!(i %in% S_ind_to_esti)){
      S_new[[counter]] = S[[i]]
      S_Vindices_new[[counter]] = S_Vindices[[i]]
      counter = counter + 1;
    }
  }
  
  res= list()
  res$S = S_new
  res$S_Vindices = S_Vindices_new
  
  return(res)
}


#' @export
get_leaves_from_node <- function(S, S_Vindices, node) {
  
  S_new = matrix(list(), 0, 1)
  S_Vindices_new = matrix(list(), 0, 1)
  
  node_indx = node$V_indices
  n = length(S)
  
  counter = 1;
  for (i in 1:n) {
    sind = S_Vindices[[i]]
    if (sum(node_indx %in% sind) > 0) {
      S_new[[counter]] = S[[i]]
      S_Vindices_new[[counter]] = S_Vindices[[i]]
      counter = counter + 1
    }
  }
  
  res= list()
  res$S = S_new
  res$S_Vindices = S_Vindices_new
  
  return(res)
}


#' @export
remove_leaves_from_node <- function(S, S_Vindices, node) {
  
  S_new = matrix(list(), 0, 1)
  S_Vindices_new = matrix(list(), 0, 1)
  
  node_indx = node$V_indices
  n = length(S)
  
  counter = 1;
  for (i in 1:n) {
    sind = S_Vindices[[i]]
    if (sum(node_indx %in% sind) == 0) {
      S_new[[counter]] = S[[i]]
      S_Vindices_new[[counter]] = S_Vindices[[i]]
      counter = counter + 1
    }
  }
  
  res= list()
  res$S = S_new
  res$S_Vindices = S_Vindices_new
  
  return(res)
}


#' @export
is_father_of_a_leaf <- function(node, S_Vindices) {
  
  is_father = FALSE
  V_ind = node$V_indices
  
  n = length(S_Vindices)
  
  for (i in 1:n) {
    q = S_Vindices[[i]]
  
    if (sum(V_ind %in% q) > 0) {
      is_father = TRUE
      return(is_father)
    }
  }
  
  return(is_father)
}


#' @export
keep_important_components <- function(S_Vindices_loop, S_Vindices, V) {
  
  n1 = length(S_Vindices_loop)
  n2 = length(S_Vindices)
  #m = length(V)
  Simp = matrix(list(), 0, 1)
  S_Vindices_imp = matrix(list(), 0, 1)
  counter = 1;
  
  for (i in 1:n1) {
    indx_loop = S_Vindices_loop[[i]]
    for (j in 1:n2) {
      indx = S_Vindices[[j]]
      if (sum(indx %in% indx_loop) > 0) {
        Simp[[counter]] = V[, indx_loop]
        S_Vindices_imp[[counter]] = indx_loop
        counter = counter + 1;
        break
      }
    }
  }
  
  res= list()
  res$S = Simp
  res$S_Vindices = S_Vindices_imp
  
  return(res)
}


#' @export
ball_line_intersection <- function(x, v, x0, R) {
  
  res = list()
  intersect = TRUE
  x = x - x0
  a = t(v) %*% v
  b = 2 * (t(x) %*% v)
  g = t(x) %*% x - R
    
  D = b^2 - 4 * a * g
    
  if (D < 0) {
      res$intersect = FALSE
      res$tmin = -Inf
      res$tmax = Inf
      return(res)
  }
  
  res$intersect = intersect
  res$tmin = (-b - sqrt(D)) / (2*a)
  res$tmax = (-b + sqrt(D)) / (2*a)
  
  return(res)
}


#' @export
count_num_of_leaves <- function(tree, num_leaves) {
  
  node = tree
  
  if (node$number_of_leaves > 0) {
    #print('HI')
    num_leaves = num_leaves + node$number_of_leaves
  }
  
  if (length(node$children) == 0) {
    return(num_leaves)
  }
  
  n = length(node$children)
  
  for (i in 1:n) {
    num_leaves = num_leaves + count_num_of_leaves(node$children[[i]], 0)
  }
  
  return(num_leaves)
}

#' @export
count_height_of_tree <- function(node) {
  
  num_leaves = 0
  while(TRUE) {
    
    print(node$number_of_leaves)
    print(length(node$children))
    if (node$number_of_leaves > 0) {
      break
    } else {
      node = node$children[[1]]
      num_leaves = num_leaves + 1
    }
  }
  
  return(num_leaves)
}


#' @export
IsInBall <- function(x, x0, r) {

  is_in = FALSE

  if (sqrt(sum((x - x0)^2)) < r) {
    is_in = TRUE
  }
  
  return(is_in)
}


#' export
find_components_new_vertices <- function(V, x0) {
  
  res = list()
  n = dim(V)[2]
  A = diag(n)
  S = matrix(list(), 0, 1)
  indices_in = c()
  cols_out = c()
  S_Vindices = matrix(list(), 0, 1)
  
  for (i in 1:n) {
    
    x = V[, i]
    if (IsInBall(x, x0, 1)) {
      cols_out = c(cols_out, i)
      next
    }
    
    for (j in (i+1):n) {
      if (j > n) {
        break
      }
      v = V[, j]
      if (IsInBall(v, x0, 1)) {
        next
      }
      v = v - x
      ball_res = ball_line_intersection(x, v, x0, 1)
      if (!ball_res$intersect) {
        A[i, j] = 1
        A[j, i] = 1
        next
      }
      if ((ball_res$tmin > 1 || ball_res$tmin < 0) && (ball_res$tmax > 1 || ball_res$tmax < 0)) {
        A[i, j] = 1
        A[j, i] = 1
      }
    }
    
  }
  
  bins = ggm::conComp(A)
  #print(bins)
  #print(cols_out)
  #print(A)
  components = unique(bins)
  #print(components)
  
  k = length(components)
  counter = 1

  for (i in 1:k) {
    if (sum(i %in% bins[cols_out]) > 0) {
      next
    }
    q = which(bins == components[i])
    S[[counter]] = V[, q]
    indices_in = c(indices_in, q)
    S_Vindices[[counter]] = q
    counter = counter + 1
  }
  
  res$S = S
  res$S_Vindices = S_Vindices
  res$V_ind_out = cols_out
  #print(cols_out)
  
  return(res)
}


#' export
check_convergence_interface <- function(A, b, V, Sind, x0, X, nu, lb, ub, last_round) {
  
  VV = V[, Sind]
  single_column = FALSE
  
  if (length(Sind) == 1){
    VV = cbind(VV, VV)
    single_column = FALSE
  }
  
  return(check_convergence(A, b, VV, x0, X, nu, lb, ub, last_round, single_column))
}


#' export
return_first_inside_interface <- function(A, b, V, cols, x0, X) {
  
  VV = V[, cols]
  single_column = FALSE
  
  if (length(cols) == 1){
    VV = cbind(VV, VV)
    single_column = FALSE
  }
  
  return(return_first_inside(A, b, VV, x0, X, single_column))
  
}


#' export
Sampling_simplex <- function(d, N) {
  
  Y = matrix( rexp(d*N, 1), d, N)
  T = colSums(Y)
  return(Y/matrix(rep(T, d), ncol=N, byrow=T))
  
}


#' export
is_in_component <- function(x, x0, A, b, V, Sind, full_check) {
  
  is_in = FALSE
  
  if (full_check) {
    q = A%*%x - b
    if (sum(q>0) > 0) {
      return(is_in)
    }
  }
  
  v = x - x0
  lambdas =  A %*% v / (b - A%*%x)
  l_max = max(lambdas)
  l_max = 1 / l_max
  #print(l_max)
  
  p = x + l_max * v
  
  Vi = V[, Sind]
  #print(Vi)
  mm = length(Sind)
  
  for (i in 1:mm) {
    
    if (mm==1){
      v = Vi - p
    } else{
      v = Vi[, i] - p
    }
    
    
    p = p - x0
    a = t(v) %*% v
    b = 2 * (t(p) %*% v)
    g = t(p) %*% p - 1
        
    D = b^2 - 4 * a * g
    #print(D)
    
    if (D < 0) {
      is_in = TRUE
      return(is_in)
    }
    
    tmin = (-b - sqrt(D)) / (2*a)
    tmax = (-b + sqrt(D)) / (2*a)
    
    if (tmin < 1 & tmin > 0) {
      next
      #is_in = false;
      #return;
    } else if(tmax < 1 & tmax > 0) {
      next
      #is_in = false;
      #return;
    } else {
      is_in = TRUE
      return(is_in)
    }
    
  }
  
  return(is_in)
}


#' export
get_fast_interior_point <- function(A, b, x0, V, Sind, cmin, cmax) {
  
  m = dim(A)[1]
  full_check = TRUE
  sc = cmax
  d = length(x0)
  
  X = boundary_randsphere(d, 2) + kronecker(matrix(1, 1, 2), matrix(x0, ncol = 1))
  x = X[, 1]
  
  while (TRUE) {
    
    bb = sc*b
    X = sample_component(A, bb, x, 1200 + floor(d^2/4), 1, x0)
    AX_b = A %*% X - kronecker(matrix(1, 1, 1200 + floor(d^2/4)), matrix(bb, ncol = 1))
    #+dim(A)[2]^2/4
    #q = colSums(sign(AX_b))
  
    #pos = which(q == -m)
    #AX_b_pos = AX_b[, pos]
  
    #pos = which(AX_b_pos == max(AX_b_pos), arr.ind = TRUE)
    pos = which(AX_b == min(colMaxs(AX_b, na.rm = TRUE)), arr.ind = TRUE)
    
    #print((A[pos[1],]%*%X[, pos[2]] - bb[pos[1]]))
    sc2 = (1.03 * (A[pos[1],]%*%X[, pos[2]] / b[pos[1]]))[1] #check this
    #print(is_in_component(X[, pos[2]], x0, A, sc2*b, V, Sind, full_check))
    x = X[, pos[2]]
    #print(sc2)
    sc=sc2
     
    if (is_in_component(x, x0, A, b, V, Sind, full_check)) {
      find_point_in_component(A, b, x0, V, cmin, cmax, x, FALSE) #just to check
      return(x)
    }
  }
  
}


#' export
get_vertices <- function(A, b) {
  
  m = dim(A)[1]
  d = m-1
  V = matrix(0, d, d+1)
  
  
  for (i in 1:m) {
    V[, i] = solve(A[-c(i),], b[-c(i)])
  }
  return(V)
}


#' export
get_distances <- function(A, b, x0) {
  
  m = dim(A)[1]
  dists = rep(0, m)
  
  for (i in 1:m) {
    dists[i] = abs(A[i,] %*% x0 - b[i]) / sqrt(sum(A[i,]^2))
  }
  #print(dists)
  #return(c(order(dists[dists>1], decreasing=TRUE), order(dists[dists<1])))
  return(order(dists, decreasing=TRUE))
}












