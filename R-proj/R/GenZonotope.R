#' Generator function for zonotopes.
#' 
#' This function can be used to generate a \eqn{d}-dimensional zonotope described by the Minkowski sum of m segments. We consider the \eqn{e_1, \dots ,e_d} generators and \eqn{m-d} random generators. Then we shift the zonotope in order to contain the origin. The origin is the center of symmetry as well. It might needs rounding before the volume computation.
#' 
#' @param dimension The dimension of the zonotope.
#' @param NumGen The number of segments that generate the zonotope.
#' 
#' @return A \eqn{m \times d} matrix that containes the \eqn{m} \eqn{d}-dimensional segments.
#' @examples 
#' # generate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' zonotope = GenZonotope(10, 20)
GenZonotope <- function(dimension, NumGen) {
  
  Zono = TRUE
  kind_gen = 0
  repr = 'zonotope'
  
  ListMat = polytope_generator(Zono, repr, kind_gen, dimension, NumGen)
  
  return(ListMat)
  
}
