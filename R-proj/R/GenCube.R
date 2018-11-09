#' Generator function for hypercubes
#' 
#' This function can be used to generate a \eqn{d}-dimensional Hypercube \eqn{[-1,1]^d} in H or V representation.
#' 
#' @param dimension The dimension of the hypercube
#' @param repr A string to declare the representation. It has to be 'H' for H-representation or 'V' for V-representation.
#' 
#' @return A hypercube in H or V-representation. For an H polytope the return value is a list with two elements: the "matrix" containing a \eqn{2d\times d} matrix \eqn{A} and the "vector" containing a \eqn{2d}-dimensional vector \eqn{b}, s.t. \eqn{Ax\leq b}. When the V-representation is chosen the return value is a \eqn{2^d \times d} matrix that containes the vertices row-wise.
#' @examples 
#' # generate a 10-dimensional hypercube in H-representation
#' PolyList = GenCube(10, 'H')
#' 
#' # generate a 15-dimension hypercube in V-representation
#' PolyList = GenCube(15, 'V')
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
  
  # remove first row
  Mat = Mat[-c(1),]
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  if (Vpoly_gen) {
    P = VPolytope(V = Mat, t=2)
  } else {
    P = HPolytope(A = -Mat, b = b, t=1)
  }
  
  return(P)
  
}
