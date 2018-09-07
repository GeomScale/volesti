#' Takes a numerical matrix in ine format and returns the matrix \eqn{A} and the vector \eqn{b} s.t.: \eqn{Ax\leq b}.
#' 
#' This function can be used to extract from a numerical matrix in ine format (see example), that describes a H-polytope, the \eqn{m\times d} matrix \eqn{A} and the \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}.
#'
#' @param A The numerical matrix in ine format (see example) of the H-polytope.
#' @return A list that containes elements "matrix" and "vector", i.e. the numerical \eqn{m\times d} matrix \eqn{A} and the numerical \eqn{m}-dimensional vector \eqn{b}, defining H-polytope \eqn{P}, s.t.:  \eqn{Ax\leq b}. For V polytopes the element "vector" is useless in practice. 
#'
#' @examples
#' # a 2d unit simplex in H-representation using numerical matrix in ine format
#' A = matrix(c(3,3,0,0,-1,0,0,0,-1,1,1,1), ncol=3, nrow=4, byrow=TRUE)
#' list_of_matrix_and_vector = modifyMat(A)
modifyMat <- function(A){
  
  # remove first row
  A = A[-c(1),]
  
  # first column is the vector b
  b = A[,1]
  
  # remove first column
  A2 = A[,-c(1)]
  
  # return final list
  retList = list("matrix"=A2, "b"=b)
  return(retList)
  
}
