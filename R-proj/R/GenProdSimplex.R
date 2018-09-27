#' Generator function for product of simplices
#' 
#' This function can be used to generate a \eqn{2d}-dimensional polytope that is defined as the product of two \eqn{d}-dimensional unit simplices in H-representation.
#' 
#' @param dimension The dimension of the simplices.
#' 
#' @return A polytope defined as the product of two unit simplices in H-representation. The retutn value is a list with two elements: the "matrix" containing a \eqn{(2d+1)\times 2d} matrix \eqn{A} and the "vector" containing a \eqn{(2d+1})-dimensional vector \eqn{b}, s.t. \eqn{Ax\leq b}.
#' @examples
#' # generate a product of two 5-dimensional simplices.
#' PolyList = GenProdSimplex(5)
#' @export
GenProdSimplex <- function(dimension) {
  
  kind_gen = 4
  m_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, m_gen)
  
  # remove first row
  Mat = Mat[-c(1),]
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = HPolytope(A = -Mat, b = b)
  
  return(P)
  
}
