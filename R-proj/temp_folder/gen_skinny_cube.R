#' Generator function for skinny hypercubes
#' 
#' This function generates a \eqn{d}-dimensional skinny hypercube \eqn{[-1,1]^{d-1}\times [-100,100]}.
#' 
#' @param dimension The dimension of the skinny hypercube.
#' 
#' @return A polytope class representing the \eqn{d}-dimensional skinny hypercube in H-representation.
#'
#' @examples
#' # generate a 10-dimensional skinny hypercube.
#' P = gen_skinny_cube(10)
#' @export
gen_skinny_cube <- function(dimension) {
  
  kind_gen = 5
  m_gen = 0
  Vpoly_gen = FALSE
  
  Mat = poly_gen(kind_gen, Vpoly_gen, FALSE, dimension, m_gen)

  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Hpolytope$new(-Mat, b, 2^(dimension -1)*200)
  
  return(P)
}
