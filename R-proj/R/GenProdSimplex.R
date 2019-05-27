#' Generator function for product of simplices
#' 
#' This function can be used to generate a \eqn{2d}-dimensional polytope that is defined as the product of two \eqn{d}-dimensional unit simplices in H-representation.
#' 
#' @param dimension The dimension of the simplices.
#' 
#' @return A polytope class representing the product of the two \eqn{d}-dimensional unit simplices in H-representation.
#' 
#' @examples
#' # generate a product of two 5-dimensional simplices.
#' P = GenProdSimplex(5)
#' @export
GenProdSimplex <- function(dimension) {
  
  kind_gen = 4
  m_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, m_gen)

  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Hpolytope$new(-Mat, b)
  
  return(P)
  
}
