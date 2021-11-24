
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
get_rad_cap <- function(x, A, b, x0) {
  
  n = x - x0
  d = 0
  
  b = b - A %*% x0
  m = dim(A)[1]
  
  row_norms = sqrt(rowSums(A^2))
  A2 = diag(1/row_norms)%*%A
  b2 = diag(1/row_norms)%*%b
  
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
  rad = sqrt(1-d^2)
  return(rad)
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
      x = y + res_int$tmin[1]*v
    } else if(res_int$tmax < 1 && res_int$tmax > 0) {
      x = y + res_int$tmax[1]*v
    } else if (res_int$tmin < 0) {
      x = y + res_int$tmin[1]*v
    } else if (res_int$tmax < 0) {
      x = y + res_int$tmax[1]*v
    } else {
      stop('false intersection')
    }
    
    if (!is_in_component(x, x0, A, b, V, S_Vindices[[i]], TRUE)) {
      stop('point outside')
    }
    
    Xs = cbind(Xs,x)
    
  }
  return(Xs)
}


#' export
get_points_on_components_imp <- function(A, b, x0, V, S, S_Vindices) {
  
  n = length(S)
  d = length(x0)
  
  y = rep(0, d)
  
  Xs = matrix(, d, 0)
  
  for (i in 1:n) {
    
    ind_verts = S_Vindices[[i]]
    n_verts = length(ind_verts)
    rad_max = 0
    p=c()
    
    for (j in 1:n_verts) {
    
      v = V[,ind_verts[j]] - y
      res_int = ball_line_intersection(y, v, x0, 1)
    
      if (!res_int$intersect) {
        stop('no intersection')
      }
      if (is_in_component(y, x0, A, b, V, ind_verts, TRUE)) {
        print("isin")
        inds = 1:n
        inds = inds[-c(i)]
        ind_verts_temp = S_Vindices[[ inds[1] ]]
        v = V[, ind_verts_temp[1]] - y
        res_int = ball_line_intersection(y, v, x0, 1)
        x = y + res_int$tmin[1]*v
      } else {
        if (res_int$tmin < 1 && res_int$tmin > 0) {
          x = y + res_int$tmin[1]*v
        } else if(res_int$tmax < 1 && res_int$tmax > 0) {
          x = y + res_int$tmax[1]*v
        } else if (res_int$tmin < 0) {
          x = y + res_int$tmin[1]*v
        } else if (res_int$tmax < 0) {
          x = y + res_int$tmax[1]*v
        } else {
          stop('false intersection')
        }
      }
    
      if (!is_in_component(x, x0, A, b, V, ind_verts, TRUE)) {
        stop('point outside')
      }
      
      rad = get_rad_cap(x, A, b, x0)
      if(rad>rad_max) {
        rad_max = rad
        p = x
      }
    }
    Xs = cbind(Xs,p)
    
  }
  return(Xs)
}


#' export
remove_small_components <- function(A, b, mus, center_2, WW, error, a_vals, ratios) {
  
  n = length(a_vals)
  log_vols = rep(0, n)
  
  a_min = Inf
  a_min_ind = 0
  
  for (i in 1:n) {
    a_seq = a_vals[[i]]
    a_min_temp = tail(a_seq, 1)
    if (a_min_temp < a_min) {
      a_min = a_min_temp
      a_min_ind = i
    }
  }
  
  inds = 1:n
  inds = inds[-c(a_min_ind)]
  ratio_min = ratios[[a_min_ind]]
  
  ratios_min = c(1)
  ratios_volumes = rep(0, n)
  ratios_volumes[a_min_ind] = 1
    
  for (i in inds) {
    
    a_seq = a_vals[[i]]
    a_max = tail(a_seq, 1)
    
    ratio_max = ratios[[i]]
    
    ratio_min_temp = prod(1/ratio_min)
    ratio_max_temp = prod(1/ratio_max)
    q = estimate_related_ratio_fischer(A, b, mus[[a_min_ind]], center_2, a_min, a_max, ratio_min_temp, ratio_max_temp, WW, error)
    
    ratios_volumes[i] = q$ratio_max_temp / (q$ratio_min_temp*ratio_min_temp)
    ratios_min = c(ratios_min, q$ratio_min_temp)
  }
  res_rem = list()
  res_rem$index = a_min_ind
  res_rem$ratios_volumes = ratios_volumes
  res_rem$ratios_min = ratios_min
  return(res_rem)
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


#' export
get_samples_2 <- function(XX, relative_vols, M) {
  
  d=dim(XX[[1]])[1]
  #print(relative_vols)
  cumsum_w = cumsum(relative_vols)
  samples = matrix(,d,0)
  for (i in 1:M) {
    r = runif(1, min = 0, max = 1)
    pos = which(cumsum_w > r)
    X = XX[[pos[1]]]
    indx <- sample(1:dim(X)[2], 1)
    samples = cbind(samples, X[,indx])
  }
  return(samples)
}



