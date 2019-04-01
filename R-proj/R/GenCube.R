#' Generator function for hypercubes
#' 
#' This function can be used to generate the \eqn{d}-dimensional unit hypercube \eqn{[-1,1]^d} in H- or V-representation.
#' 
#' @param dimension The dimension of the hypercube
#' @param repr A string to declare the representation. It has to be \code{'H'} for H-representation or \code{'V'} for V-representation.
#' 
#' @return A polytope class representing the unit \eqn{d}-dimensional hypercube in H- or V-representation.
#' @examples 
#' # generate a 10-dimensional hypercube in H-representation
#' P = GenCube(10, 'H')
#' 
#' # generate a 15-dimension hypercube in V-representation
#' P = GenCube(15, 'V')
#' @export
GenCube <- function(dimension, repr) {
  
  kind_gen = 1
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
