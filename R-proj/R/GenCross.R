#' Generator function for cross polytopes.
#' 
#' This function can be used to generate a d-dimensional cross polytope in H or V representation.
#' 
#' @param Dimension The dimension of the cross polytope.
#' @param represenation A string to declare the representation. It has to be H for H-representation or V for V-representation.
#' 
#' @return A cross polytope in H or V-representation. For an H polytope the return value is a list with two elements: the "matrix" containing a \eqn{2^d \times d} matrix A and the "vector" containing a \eqn{2^d} -dimensional vector b, s.t. \eqn{Ax\leq b}. When the V-representation is chosen the return value is a \eqn{2d \times d} matrix that containes the vertices row-wise.
#' @examples 
#' # generate a 10-dimensional cross polytope in H-representation
#' PolyList = GenCross(10, 'H')
#' 
#' # generate a 15-dimension cross polytope in V-representation
#' PolyList = GenCross(15, 'V')
GenCross <- function(DimGen, repr) {
  
  Zono = FALSE
  kind_gen = 2
  NumGen = 0
  
  ListMat = polytope_generator(Zono, repr, kind_gen, DimGen, NumGen)
  
  return(ListMat)
  
}