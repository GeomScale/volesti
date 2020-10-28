#' Generator function for random V-polytopes
#' 
#' This function generates a \eqn{d}-dimensional polytope in V-representation with \eqn{m} vertices. We pick \eqn{m} random points from the boundary of the \eqn{d}-dimensional unit hypersphere as vertices.
#' 
#' @param dimension The dimension of the convex polytope.
#' @param nvertices The number of the vertices.
#' @param generator The body that the generator samples uniformly the vertices from: (a) 'cube' or (b) 'sphere'.
#' @param seed Optional. A fixed seed for the generator.
#' 
#' @return A polytope class representing a V-polytope.
#' @examples 
#' # generate a 10-dimensional polytope defined as the convex hull of 25 random vertices
#' P = gen_rand_vpoly(10, 25)
#' @export
gen_rand_vpoly <- function(dimension, nvertices, generator = NULL, seed = NULL) {
  
  kind_gen = 4
  
  if(!is.null(generator)){
    if (generator == 'cube'){
      kind_gen = 5
    } else if (generator != 'sphere') {
      stop("Wrong generator!")
    }
  }
  
  Mat = poly_gen(kind_gen, TRUE, FALSE, dimension, nvertices, seed)

  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Vpolytope$new(Mat)
  
  return(P)
}
