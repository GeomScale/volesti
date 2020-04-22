#' Generator function for skinny hypercubes
#' 
#' This function can be used to generate a \eqn{d}-dimensional skinny hypercube in H-representation.
#' 
#' @param dimension The dimension of the skinny hypercube.
#' 
#' @return A polytope class representing the \eqn{d}-dimensional skinny hypercube in H-representation.
#'
#' @examples
#' # generate a 10-dimensional skinny hypercube.
#' P = GenSkinnyCube(10)
#' @export
GenSkinnyCube <- function(dimension) {
  
  kind_gen = 5
  m_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, m_gen)

  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Hpolytope$new(-Mat, b)
  
  return(P)
}
