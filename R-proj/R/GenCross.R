#' Generator function for cross polytopes
#' 
#' This function can be used to generate the \eqn{d}-dimensional cross polytope in H- or V-representation.
#' 
#' @param dimension The dimension of the cross polytope.
#' @param repr A string to declare the representation. It has to be \code{'H'} for H-representation or \code{'V'} for V-representation.
#' 
#' @return A polytope class representing a cross polytope in H- or V-representation.
#' @examples 
#' # generate a 10-dimensional cross polytope in H-representation
#' P = GenCross(10, 'H')
#' 
#' # generate a 15-dimension cross polytope in V-representation
#' P = GenCross(15, 'V')
#' @export
GenCross <- function(dimension, repr) {
  
  kind_gen = 2
  m_gen = 0
  if (repr == "V") {
    Vpoly_gen = TRUE
  } else if (repr == "H") {
    Vpoly_gen = FALSE
  } else {
    stop('Not a known representation.')
  }
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, m_gen)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  if (Vpoly_gen) {
    P = Vpolytope$new(Mat)
  } else {
    P = Hpolytope$new(-Mat, b)
  }
  
  return(P)
}