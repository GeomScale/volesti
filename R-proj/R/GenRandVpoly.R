#' Generator function for random V-polytopes
#' 
#' This function can be used to generate a \eqn{d}-dimensional polytope in V-representation with \eqn{m} vertices. We pick \eqn{m} random points from the boundary of the \eqn{d}-dimensional unit hypersphere as vertices.
#' 
#' @param dimension The dimension of the convex polytope.
#' @param m The number of the vertices.
#' 
#' @return A polytope class representing a V-polytope.
#' @examples 
#' # generate a 10-dimensional polytope defined as the convex hull of 25 random vertices
#' P = GenRandVpoly(10, 25)
#' @export
GenRandVpoly <- function(dimension, m, body = NULL) {
  
  if (is.null(body) || body == 'sphere'){
    kind_gen = 4
  } else if (body == 'cube') {
    kind_gen = 5
  } else {
    warning("Generators use only random points from the unit cube or the unit sphere")
    return(Vpolytope$new())
  }
  Vpoly_gen = TRUE
  
  Mat = poly_gen(kind_gen, FALSE, Vpoly_gen, dimension, m)

  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Vpolytope$new(Mat)
  
  return(P)
}
