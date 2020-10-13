#' Generator function for random H-polytopes
#' 
#' This function generates a \eqn{d}-dimensional polytope in H-representation with \eqn{m} facets. We pick \eqn{m} random hyperplanes tangent on the \eqn{d}-dimensional unit hypersphere as facets.
#' 
#' @param dimension The dimension of the convex polytope.
#' @param nfacets The number of the facets.
#' @param generator A list that could contain two elements: (a) 'generator', the constants \eqn{b_i} for each facets from: (i) 'sphere' from the surface of the unit hypersphere, (ii) 'ball' the interior of the unit hypersphere, the defalut value is 'sphere', and (b) 'seed' to set a spesific seed for the number generator.
#' 
#' @return A polytope class representing a H-polytope.
#' @examples 
#' # generate a 10-dimensional polytope with 50 facets
#' P = gen_rand_hpoly(10, 50)
#' @export
gen_rand_hpoly <- function(dimension, nfacets, generator = list('generator' = 'sphere')) {
  
  seed = NULL
  if (!is.null(generator$seed)) {
    seed = generator$seed
  }
  
  if (is.null(generator$generator)) {
    kind_gen = 6
  } else if (generator$generator == 'ball'){
    kind_gen = 7
  } else if (generator$generator == 'sphere') {
    kind_gen = 6
  } else {
    stop("Wrong generator!")
  }
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, FALSE, dimension, nfacets, seed)
  
  # first column is the vector b
  b = Mat[, 1]
  Mat = Mat[, c(1), drop = FALSE]
  
  P = Hpolytope(A = Mat, b = b)
  
  return(P)
}
