
#' export
get_max_cap <- function(X, A, b, x0, V, Sind, N){
  
  N0 = N
  indxs = c()
  while (TRUE) {
    indxs = unique(c(indxs, sample(1:dim(X)[2], N0)))
    if (length(indxs) == N) {
      break
    } else {
      N0 = N - length(indxs)
    }
  }
  
  rad = 2
  best_point = c()
  b = b - A %*% x0
  m = dim(A)[1]
  
  row_norms = sqrt(rowSums(A^2))
  A2 = diag(1/row_norms)%*%A
  b2 = diag(1/row_norms)%*%b
  
  for (i in 1:N) {
    
    n = X[, indxs[i]] - x0
    d = 0
    
    for (j in 1:m) {
      
      bi = b2[j]
      ai = A2[j,]
      
      ai_n = ai%*%n
      ai_n = ai_n[1]
      
      a = -ai_n*ai_n - 1
      bb = 2 * bi * ai_n
      g = 1 - bi*bi
      
      D = bb^2 - 4*a*g
      
      if (D<0) {
        next
      }
      
      d1 = (-bb - sqrt(D)) / (2*a)
      d2 = (-bb + sqrt(D)) / (2*a)
      
      if (d1 > d) {
        d = d1
      }
      if (d2 > d) {
        d = d1
      }
      
    }
    
    if (d < rad) {
      rad = d
      best_point = n + x0
    }
    
  }
  
  rad = sqrt(1-rad^2)
  res = list()
  res$rad = rad
  res$center = best_point
  return(res)
  
}


#' export
get_points_on_components <- function(A, b, x0, V, S, S_Vindices, V_ind_out) {
  
  n = length(S)
  d = length(x0)
  
  y = rep(0, d)
  
  Xs = matrix(, d, 0)
  
  for (i in 1:n) {
    
    Verts = V[,S_Vindices[[i]]]
    if (is.null(dim(Verts))) {
      v = Verts - y
    } else {
      v = Verts[, 1] - y
    }
    res_int = ball_line_intersection(y, v, x0, 1)
    
    if (!res_int$intersect) {
      stop('no intersection')
    }
    
    if (res_int$tmin < 1 && res_int$tmin > 0) {
      x = y + res_int$tmin*v
    } else if(res_int$tmax < 1 && res_int$tmax > 0) {
      x = y + res_int$tmax*v
    } else if (res_int$tmin < 0) {
      x = y + res_int$tmin*v
    } else if (res_int$tmax < 0) {
      x = y + res_int$tmax*v
    } else {
      stop('false intersection')
    }
    
    if (!is_in_component(x, x0, A, b, V, S_Vindices[[i]], TRUE)) {
      stop('point outside')
    }
    
    Xs = cbind(Xs,x)
    
  }
  return(XS)
}


#' export
remove_small_components <- function(a_vals, ratios, d) {
  
  n = length(a_vals)
  log_vols = rep(0, n)
  
  for (i in 1:n) {
    
    a_seq = a_vals[[i]]
    a_max = tail(a_seq, 1)
    rats = ratios[[i]]
    
    log_vols[i] = ((d/2) * log(2*pi) + log(compute_bessel(d, a_max)) - (d/2 - 1)*a_max) - sum(log(rats))
    
  }
  
}


#' export
log_integral_fischer <-function(d, k) {
  
  a = compute_bessel(d/2-1, k)
  if (is.nan(a) || is.infinite(a)) {
    return(NaN)
  }
  log_int = (d/2) * log(2*pi) + log(a) - (d/2-1)*log(k)
  return(log_int)
}


#' export
log_gamma_function <- function(x) {
  
  if (x<=100){
    return(log(gamma(x)))
  }
  
  return(log(x-1) + log_gamma_function(x-1))
  
}


#' export
log_volume_n_sphere <- function(d) {
  
  log_vol = (d/2)*log(2*pi) - log_gamma_function(d/2)
  
  return(log_vol)
}




