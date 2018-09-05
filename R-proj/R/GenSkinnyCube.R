#' Generator function for skinny hypercubes.
#' 
#' This function can be used to generate a d-dimensional skinny hypercube in H-representation.
#' 
#' @param dimension The dimension of the skinny hypercube.
#' 
#' @return A d-dimensional skinny hypercube in H-representation.The retutn value is a list with two elements: the "matrix" containing a \eqn{2d\times d} matrix A and the "vector" containing a 2d-dimensional vector b, s.t. \eqn{Ax\leq b}.
#' @examples
#' # generate a 10-dimensional skinny hypercube.
#' PolyList = GenSkinnyCube(10)
GenSkinnyCube <- function(dimension, repr = 'H') {
  
  Zono = FALSE
  kind_gen = 5
  NumGen = 0
  
  ListMat = polytope_generator(Zono, repr, kind_gen, dimension, NumGen)
  
  return(ListMat)
  
}
