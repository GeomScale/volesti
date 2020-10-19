#' Generator function for random V-polytopes
#' 
#' This function generates a \eqn{d}-dimensional polytope in V-representation with \eqn{m} vertices. We pick \eqn{m} random points from the boundary of the \eqn{d}-dimensional unit hypersphere as vertices.
#' 
#' @param dimension The dimension of the convex polytope.
#' @param nvertices The number of the vertices.
#' @param generator A list that could contain two elements.
#' \itemize{
#' \item{body }{ the body that the generator samples uniformly the vertices from: (i) 'cube' or (ii) 'sphere', the default value is 'sphere'.}
#' \item{seed }{ Optional. A fixed seed for the number generator.}
#' }
#' 
#' @return A polytope class representing a V-polytope.
#' @examples 
#' # generate a 10-dimensional polytope defined as the convex hull of 25 random vertices
#' P = gen_rand_vpoly(10, 25)
#' @export
gen_rand_vpoly <- function(dimension, nvertices, generator = list('body' = 'sphere')) {
  
  seed = NULL
  if (!is.null(generator$seed)) {
    seed = generator$seed
  }
  
  if (is.null(generator$body)) {
    kind_gen = 4
  } else if (generator$body == 'cube'){
    kind_gen = 5
  } else if (generator$body == 'sphere') {
    kind_gen = 4
  } else {
    stop("Wrong generator!")
  }
  
  Mat = poly_gen(kind_gen, TRUE, FALSE, dimension, nvertices, seed)

  # first column is the vector b
  b = Mat[, 1]
  Mat = Mat[, -c(1), drop = FALSE]
  
  P = Vpolytope(V = Mat)
  
  return(P)
}
