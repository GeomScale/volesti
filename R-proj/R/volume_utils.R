
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