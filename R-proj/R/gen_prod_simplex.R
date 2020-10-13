#' Generator function for product of simplices
#' 
#' This function generates a \eqn{2d}-dimensional polytope that is defined as the product of two \eqn{d}-dimensional unit simplices in H-representation.
#' 
#' @param dimension The dimension of the simplices.
#' 
#' @return A polytope class representing the product of the two \eqn{d}-dimensional unit simplices in H-representation.
#' 
#' @examples
#' # generate a product of two 5-dimensional simplices.
#' P = gen_prod_simplex(5)
#' @export
gen_prod_simplex <- function(dimension) {
  
  kind_gen = 4
  m_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, FALSE, dimension, m_gen)

  # first column is the vector b
  b = Mat[, 1]
  Mat = Mat[, -c(1), drop = FALSE]
  
  P = Hpolytope(A = -Mat, b = b, volume = (1/prod(1:dimension))^2)
  
  return(P)
  
}
