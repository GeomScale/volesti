#' Generator function for skinny hypercubes
#' 
#' This function can be used to generate a \eqn{d}-dimensional skinny hypercube only in H-representation.
#' 
#' @param dimension The dimension of the skinny hypercube.
#' 
#' @return A \eqn{d}-dimensional skinny hypercube in H-representation. The retutn value is a list with two elements: the "matrix" containing a \eqn{2d\times d} matrix \eqn{A} and the "vector" containing a \eqn{2d}-dimensional vector \eqn{b}, s.t. \eqn{Ax\leq b}.
#' @examples
#' # generate a 10-dimensional skinny hypercube.
#' PolyList = GenSkinnyCube(10)
#' @export
GenSkinnyCube <- function(dimension) {
  
  kind_gen = 5
  m_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, m_gen)
  
  # remove first row
  Mat = Mat[-c(1),]
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = HPolytope$new(-Mat, b)
  
  return(P)
}
