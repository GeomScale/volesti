
#' export
compute_cap_volume <- function(A, b, x0, x, V, indx) {
  
  d = length(x0)
  res = get_hyperplane(V, indx, x0)
  n = res$n
  bi = res$bi
  
  rad = sqrt(1 - bi^2)
  cen = bi*n
  
  #y1 = sample_direction_2(n, rad) + cen
  y2 = sample_direction_2(n, rad) + cen
  
  a = acos((n%*%y2)[1])
  
  N = 200000
  X = runif(N, min = 0, max = a)
  
  vol = exp(log_volume_n_sphere(d-1))
  sc = sin(X)^(d-2)
  vol = vol * (mean(sc)*a)
  return(vol)
}


#' export
sample_direction_2 <- function(n, rad) {
  
  d = length(n)
  ut = rnorm(d, mean=0, sd=1)
  ut = ut / sqrt(sum(ut^2))
  
  u = (diag(d) - n%*%t(n))%*%ut# / norm((diag(n) - n%*%t(n))%*%ut)
  u = (rad*u) / sqrt(sum(u^2))
  #u = rad * u
  return(u)
}


#' export
get_hyperplane <- function(V, indx, x0) {
  
  d = length(x0)
  indices = 1:(d+1)
  ind_out = indices[-c(indx)]
  
  x = V[,indx]
  Y = matrix(, 0, d)
  for (i in 1:length(ind_out)) {
    
    v = V[,ind_out[i]] - x
    ball_res = ball_line_intersection(x, v, x0, 1)
    
    y = x + ball_res$tmin[1]*v
    Y = rbind(Y, t(y))
    
  }
  
  b = rep(1, d)
  n = solve(Y, b)
  
  b = b / sqrt(sum(n^2))
  n = n / sqrt(sum(n^2))
  
  b = b - (n%*%x0)[1]
  
  res = list()
  res$n = n
  res$bi = b[1]
  
  return(res)
}



