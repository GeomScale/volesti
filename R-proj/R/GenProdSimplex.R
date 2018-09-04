#' Generator function for product of simplices.
#' 
#' This function can be used to generate a 2d-dimensional polytope that is defined as the product of two d-dimensional unit simplices.
#' 
#' @param Dimension The dimension of the simplex.
#' 
#' @return A polytope defined as the product of two unit simplices in H-representation.The retutn value is a list with two elements: the "matrix" containing a \eqn{2d+1\times 2d} matrix A and the "vector" containing a 2d+1-dimensional vector b, s.t. \eqn{Ax\leq b}.
#' @examples
#' # generate a product of two 5-dimensional simplices.
#' PolyList = GenProdSimplex(5)
GenProdSimplex <- function(DimGen, repr = 'H') {
  
  Zono = FALSE
  kind_gen = 4
  NumGen = 0
  
  ListMat = polytope_generator(Zono, repr, kind_gen, DimGen, NumGen)
  
  return(ListMat)
  
}
