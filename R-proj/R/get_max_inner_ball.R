get_max_inner_ball <- function(A, b, sparseness = FALSE) {
  
  m = dim(A)[1]
  d = dim(A)[2]
  
  row_ind = c()
  col_ind = c()
  values = c()
  
  if (sparseness){
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
  } else {
    for (i in 1:m) {
      row_ind = c(row_ind, rep(i, d+1))
    }
    col_ind = rep(1:(d+1), m)
    A = cbind(A, sqrt(rowSums(A^2)))
    values = as.vector(t(A))
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
  
  #print(paste0("sol = ", r$sol$itr$xx))
  
  ret_list = list()
  ret_list$center = r$sol$itr$xx[1:d]
  ret_list$radius = radius
  
  if (radius <=0 ){
    #stop("negative radius")
  }
  
  return(ret_list)
  
}
  