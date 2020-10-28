#' Generator function for zonotopes
#' 
#' This function generates a random \eqn{d}-dimensional zonotope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
#' The function considers \eqn{m} random directions in \eqn{R^d}. There are three strategies to pick the length of each segment: a) it is uniformly sampled from \eqn{[0,100]}, b) it is random from \eqn{\mathcal{N}(50,(50/3)^2)} truncated to \eqn{[0,100]}, c) it is random from \eqn{Exp(1/30)} truncated to \eqn{[0,100]}.
#' 
#' @param dimension The dimension of the zonotope.
#' @param nsegments The number of segments that generate the zonotope.
#' @param generator The distribution to pick the length of each segment from \eqn{[0,100]}: (a) 'uniform', (b) 'gaussian' or (c) 'exponential'.
#' @param seed Optional. A fixed seed for the generator.
#'  
#' @return A polytope class representing a zonotope.
#'
#' @examples 
#' # generate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' P = gen_rand_zonotope(10, 20)
#' @export
gen_rand_zonotope <- function(dimension, nsegments, generator = NULL, seed = NULL) {
  
  kind_gen = 1
  
  if (!is.null(generator)) {
    if (generator == 'gaussian') {
      kind_gen = 2
    } else if (generator == 'exponential') {
      kind_gen = 3
    } else if (generator != 'uniform'){
      stop("Wrong generator!")
    }
  }
  
  Mat = poly_gen(kind_gen, FALSE, TRUE, dimension, nsegments, seed)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Zonotope$new(Mat)

  return(P)
}
