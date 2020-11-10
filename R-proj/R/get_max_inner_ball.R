get_max_inner_ball <- function(A, b) {
  
  m = dim(A)[1]
  d = dim(A)[2]
  
  row_ind = c()
  col_ind = c()
  values = c()
  
  for (i in 1:m) {
    for (j in 1:d) {
      if (A[i, j] != 0) {
        row_ind = c(row_ind, i)
        col_ind = c(col_ind, j)
        values = c(values, A[i, j])
      }
    }
    row_ind = c(row_ind, i)
    col_ind = c(col_ind, d+1)
    values = c(values, sqrt(sum(A[i,]^2)))
  }
  
  prob <- list()
  prob$A <- sparseMatrix(row_ind, col_ind, x=values)
  
  # Bound values for constraints
  prob$bc <- rbind(blc=c(rep(-Inf, m)), 
                   buc=c(b))
  
  # Bound values for variables
  prob$bx <- rbind(blx = rep(-Inf, d+1), 
                   bux = rep(Inf, d+1))
  
  prob$sense <- "max"
  obj = rep(0, d+1)
  obj[d+1] = 1
  prob$c = obj
  r <- mosek(prob, list(verbose=0))
  stopifnot(identical(r$response$code, 0))
  radius = prob$c %*% r$sol$itr$xx
  
  print(paste0("sol = ", r$sol$itr$xx))
  
  ret_list = list()
  ret_list$center = r$sol$itr$xx[1:d]
  ret_list$radius = radius
  
  return(ret_list)
  
}
  