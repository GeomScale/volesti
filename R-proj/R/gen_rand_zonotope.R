#' Generator function for zonotopes
#' 
#' This function can be used to generate a random \eqn{d}-dimensional zonotope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments. We consider \eqn{m} random directions in \eqn{R^d} and for each direction we pick a random length in \eqn{[(,\sqrt{d}]} in order to define \eqn{m} segments.
#' 
#' @param dimension The dimension of the zonotope.
#' @param n_segments The number of segments that generate the zonotope.
#' @param generator The distribution to pick the length of each segment from \eqn{[0,100]}: (a) 'uniform', (b) 'gaussian' or (c) 'exponential'.
#' 
#' @return A polytope class representing a zonotope.
#'
#' @examples 
#' # generate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' P = gen_rand_zonotope(10, 20)
#' @export
gen_rand_zonotope <- function(dimension, n_segments, generator = NULL) {
  
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
  
  Mat = poly_gen(kind_gen, FALSE, TRUE, dimension, n_segments)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Zonotope$new(Mat)

  return(P)
}
