#' Generator function for hypercubes.
#' 
#' This function can be used to generate a d-dimensional Hypercube \eqn{[-1,1]^d} in H or V representation.
#' 
#' @param dimension The dimension of the hypercube
#' @param repr A string to declare the representation. It has to be H for H-representation or V for V-representation.
#' 
#' @return A hypercube in H or V-representation. For an H polytope the return value is a list with two elements: the "matrix" containing a \eqn{2d\times d} matrix A and the "vector" containing a 2d-dimensional vector b, s.t. \eqn{Ax\leq b}. When the V-representation is chosen the return value is a \eqn{2^d \times d} matrix that containes the vertices row-wise.
#' @examples 
#' # generate a 10-dimensional hypercube in H-representation
#' PolyList = GenCube(10, 'H')
#' 
#' # generate a 15-dimension hypercube in V-representation
#' PolyList = GenCube(15, 'V')
GenCube <- function(dimension, repr) {
  
  Zono = FALSE
  kind_gen = 1
  NumGen = 0
  
  ListMat = polytope_generator(Zono, repr, kind_gen, dimension, NumGen)
  
  return(ListMat)
  
}
