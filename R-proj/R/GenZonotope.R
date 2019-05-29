#' Generator function for zonotopes
#' 
#' This function can be used to generate a random \eqn{d}-dimensional zonotope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments. We consider \eqn{m} random directions in \eqn{R^d} and for each direction we pick a random length in \eqn{[(,\sqrt{d}]} in order to define \eqn{m} segments.
#' 
#' @param dimension The dimension of the zonotope.
#' @param NumGen The number of segments that generate the zonotope.
#' 
#' @return A polytope class representing a zonotope.
#'
#' @examples 
#' # generate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' P = GenZonotope(10, 20)
#' @export
GenZonotope <- function(dimension, NumGen) {
  
  kind_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, NumGen)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Zonotope$new(Mat)

  return(P)
}
