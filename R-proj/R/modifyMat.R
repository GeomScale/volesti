#' Takes a numerical matrix in ine format and returns the matrix A and the vector b.
#' 
#' This function can be used to extract from a numerical matrix in ine format (see example), that describes a H-polytope, the \eqn{m\times d} matrix A and the d-dimensional vector b, s.t.: \eqn{Ax\leq b}.
#'
#' @param A The numerical matrix in ine format (see example) of the H-polytope.
#' @return A list that contains the numerical \eqn{m\times d} matrix A and the numerical d-dimensional vector b, defining H-polytope P, s.t.:  \eqn{Ax\leq b}.
#'
#' @examples
#' # a 2d unit simplex in H-representation using numerical matrix in ine format
#' A = matrix(c(3,3,0,0,-1,0,0,0,-1,1,1,1), ncol=3, nrow=4, byrow=TRUE)
#' list_of_matrix_and_vector = modifyMat(A)
modifyMat <- function(A){
  
  b=A[,1]
  A2=A[,-c(1)]
  retList=list("matrix"=A2, "vector"=b)
  return(retList)
  
}

