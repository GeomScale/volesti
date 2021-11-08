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
    if (!ball_line_res$intersection || ((ball_line_res$tmin>1 || ball_line_res$tmin<0) & (ball_line_res$tmax>1 || ball_line_res$tmax<0))) {
  
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
  node$isLeaf = false
  node$number_of_leaves = 0
  node$ratio_of_leaves = c()
  node$S_ind_to_esti = c()
  node$S_Vindices_leaves = matrix(list(), 0, 1)
  node$children = matrix(list(), 0, 1)
  
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
      S_new{counter} = S[[i]]
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
      S_new{counter} = S[[i]]
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
  m = length(V)
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
  
  if (node.number_of_leaves > 0) {
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
      
      v = V[, j]
      if (IsInBall(v, x0, 1)) {
        next
      }
      v = v - x
      ball_res = ball_line_intersection(x, v, x0, 1)
      if (!ball_res$intersection) {
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
  components = unique(bins)
  
  k = length(components)
  counter = 1

  for (i in 1:k) {
    if (sum(i %in% bins[cols_out]) > 0) {
      next
    }
    q = which(bins == components(i))
    S[[counter]] = V[, q]
    indices_in = c(indices_in, q)
    S_Vindices[[counter]] = q[q<=n]
    counter = counter + 1
  }
  
  res$S = S
  res$S_Vindices = S_Vindices
  res$V_ind_out = cols_out
  
  return(res)
}

