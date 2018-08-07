#' takes a numerical matrix in ine format and return numerical matrix A and vector b: Ax<=b
#'
#' @param A the numerical matrix in ine format of the H-polytope
#' @return numerical matrix A and vector b: Ax<=b
#'
#' @examples
#' modifyMat(A)
modifyMat <- function(A){
  
  b=A[,1]
  A2=A[,-c(1)]
  retList=list("matrix"=A2, "vector"=b)
  return(retList)
  
}

