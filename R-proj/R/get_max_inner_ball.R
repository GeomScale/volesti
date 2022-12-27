get_max_inner_ball <- function(A, b) {
  
  d = dim(A)[2]

  A = -t(cbind(A, sqrt(rowSums(A^2))))
  b = -b

  dvec = t(rep(0,d+1))
  dvec[d+1] = 1
  
  Dmat = diag(d+1) * 1e-10

  ret = quadprog::solve.QP(Dmat, dvec, A, b)
  
  ret_list = list()
  ret_list$center = ret$solution[1:d]
  ret_list$radius = ret$solution[d+1]
  
  if (ret_list$radius <=0 ){
    stop("negative radius")
  }
  
  return(ret_list)
  
}
