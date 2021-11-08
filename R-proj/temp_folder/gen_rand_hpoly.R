#' Generator function for random H-polytopes
#' 
#' This function generates a \eqn{d}-dimensional polytope in H-representation with \eqn{m} facets. We pick \eqn{m} random hyperplanes tangent on the \eqn{d}-dimensional unit hypersphere as facets.
#' 
#' @param dimension The dimension of the convex polytope.
#' @param nfacets The number of the facets.
#' @param seed Optional. A fixed seed for the generator.
#' 
#' @return A polytope class representing a H-polytope.
#' @examples 
#' # generate a 10-dimensional polytope with 50 facets
#' P = gen_rand_hpoly(10, 50)
#' @export
gen_rand_hpoly <- function(dimension, nfacets, seed = NULL) {
  
  kind_gen = 6
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, FALSE, dimension, nfacets, seed)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Hpolytope$new(Mat, b)
  
  return(P)
}
