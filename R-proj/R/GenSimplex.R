#' Generator function for simplices
#' 
#' This function can be used to generate a \eqn{d}-dimensional unit simplex in H or V representation.
#' 
#' @param dimension The dimension of the simplex.
#' @param repr A string to declare the representation. It has to be 'H' for H-representation or 'V' for V-representation.
#' 
#' @return A simplex in H or V-representation. For an H polytope the return value is a list with two elements: the "matrix" containing a \eqn{(d+1)\times d} matrix \eqn{A} and the "vector" containing a \eqn{(d+1)}-dimensional vector \eqn{b}, s.t. \eqn{Ax\leq b}. When the V-representation is chosen the return value is a \eqn{(d+1)\times d} matrix that containes the vertices row-wise.
#' @examples
#' # generate a 10-dimensional simplex in H-representation
#' PolyList = GenSimplex(10, 'H')
#' 
#' # generate a 20-dimensional simplex in V-representation
#' PolyList = GenSimplex(20, 'V')
#' @export
GenSimplex <- function(dimension, repr) {
  
  kind_gen = 3
  m_gen = 0
  if (repr == "V") {
    Vpoly_gen = TRUE
  } else if (repr == "H") {
    Vpoly_gen = FALSE
  } else {
    stop('Not a known representation.')
  }
  
  Mat = poly_gen(kind_gen, Vpoly_gen, dimension, m_gen)
  
  # remove first row
  Mat = Mat[-c(1),]
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  if (Vpoly_gen) {
    P = VPolytope$new(Mat)
  } else {
    P = HPolytope$new(-Mat, b)
  }
  
  return(P)
}
