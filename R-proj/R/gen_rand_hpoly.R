#' Generator function for random H-polytopes
#' 
#' This function generates a \eqn{d}-dimensional polytope in H-representation with \eqn{m} facets. We pick \eqn{m} random hyperplanes tangent on the \eqn{d}-dimensional unit hypersphere as facets.
#' 
#' @param dimension The dimension of the convex polytope.
#' @param nfacets The number of the facets.
#' @param generator A list that could contain two elements.
#' \itemize{
#' \item{constants }{ To declare how to set the constants \eqn{b_i} for each facets: (i) 'sphere', each hyperplane is tangent to the hypersphere of radius 10, (ii) 'uniform' for each \eqn{b_i} the generator picks a uniform number from \eqn{(0,1)}. The defalut value is 'sphere'.}
#' \item{seed }{ Optional. A fixed seed for the number generator.}
#' }
#'  
#' @return A polytope class representing a H-polytope.
#' @examples 
#' # generate a 10-dimensional polytope with 50 facets
#' P = gen_rand_hpoly(10, 50)
#' @export
gen_rand_hpoly <- function(dimension, nfacets, generator = list('constants' = 'sphere')) {
  
  seed = NULL
  if (!is.null(generator$seed)) {
    seed = generator$seed
  }
  
  if (is.null(generator$constants)) {
    kind_gen = 6
  } else if (generator$constants == 'sphere'){
    kind_gen = 6
  } else if (generator$constants == 'uniform') {
    kind_gen = 7
  } else {
    stop("Wrong generator!")
  }
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, FALSE, dimension, nfacets, seed)
  
  # first column is the vector b
  b = Mat[, 1]
  Mat = Mat[, -c(1), drop = FALSE]
  
  P = Hpolytope(A = Mat, b = b)
  
  return(P)
}
