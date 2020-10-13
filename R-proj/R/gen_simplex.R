#' Generator function for simplices
#' 
#' This function generates the \eqn{d}-dimensional unit simplex in H- or V-representation.
#' 
#' @param dimension The dimension of the unit simplex.
#' @param representation A string to declare the representation. It has to be \code{'H'} for H-representation or \code{'V'} for V-representation.  Default valus is 'H'.
#' 
#' @return A polytope class representing the \eqn{d}-dimensional unit simplex in H- or V-representation.
#' @examples
#' # generate a 10-dimensional simplex in H-representation
#' PolyList = gen_simplex(10, 'H')
#' 
#' # generate a 20-dimensional simplex in V-representation
#' P = gen_simplex(20, 'V')
#' @export
gen_simplex <- function(dimension, representation = 'H') {
  
  kind_gen = 3
  m_gen = 0
  if (representation == "V") {
    Vpoly_gen = TRUE
  } else if (representation == "H") {
    Vpoly_gen = FALSE
  } else {
    stop('Not a known representation.')
  }
  
  Mat = poly_gen(kind_gen, Vpoly_gen, FALSE, dimension, m_gen)

  # first column is the vector b
  b = Mat[, 1]
  Mat = Mat[, -c(1), drop = FALSE]
  
  if (Vpoly_gen) {
    P = Vpolytope(V = Mat, volume = 1/prod(1:dimension))
  } else {
    P = Hpolytope(A = -Mat, b = b, volume = 1/prod(1:dimension))
  }
  
  return(P)
}
