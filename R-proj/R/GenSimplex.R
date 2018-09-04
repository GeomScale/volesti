#' Generator function for simplices.
#' 
#' This function can be used to generate a d-dimensional unit simplex in H or V representation.
#' 
#' @param Dimension The dimension of the simplex.
#' @param representation A string to declare the representation. It has to be H for H-representation or V for V-representation.
#' 
#' @return A simplex in H or V-representation. For an H polytope the return value is a list with two elements: the "matrix" containing a \eqn{d+1\times d} matrix A and the "vector" containing a d+1-dimensional vector b, s.t. \eqn{Ax\leq b}. When the V-representation is chosen the return value is a \eqn{d+1\times d} matrix that containes the vertices row-wise.
#' @examples
#' # generate a 10-dimensional simplex in H-representation
#' PolyList = GenSimplex(10, 'H')
#' 
#' # generate a 20-dimensional simplex in V-representation
#' PolyList = GenSimplex(20, 'V')
GenSimplex <- function(DimGen, repr) {
  
  Zono = FALSE
  kind_gen = 3
  NumGen = 0
  
  ListMat = polytope_generator(Zono, repr, kind_gen, DimGen, NumGen)
  
  return(ListMat)
  
}
