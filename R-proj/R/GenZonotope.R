#' Generator function for zonotopes
#' 
#' This function can be used to generate a \eqn{d}-dimensional zonotope defined by the Minkowski sum of \eqn{m} segments. We consider the \eqn{e_1, \dots ,e_d} generators and \eqn{m-d} random generators. Then we shift the zonotope in order to contain the origin. The origin is the center of symmetry as well. It might needs rounding before the volume approximation using SequenceOfBalls or CoolingGaussian algorithms.
#' 
#' @param dimension The dimension of the zonotope.
#' @param NumGen The number of segments that generate the zonotope.
#' 
#' @return A zonotope.
#' @examples 
#' # generate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' zonotope = GenZonotope(10, 20)
#' @export
GenZonotope <- function(dimension, NumGen) {
  
  kind_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, NumGen)
  
  # remove first row
  Mat = Mat[-c(1),]
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Zonotope$new(Mat)

  return(P)
}
