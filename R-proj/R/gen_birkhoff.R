#' Generator function for Birkhoff polytope
#' 
#' This function can be used to generate the full dimensional \eqn{n}-Birkhoff polytope in H-representation.
#' The dimension of the generated polytope is \eqn{(n-1)^2}.
#' 
#' @param n The order of the Birkhoff polytope
#' 
#' @return A polytope class representing the full dimensional \eqn{n}-Birkhoff polytope in H-representation.
#' @examples 
#' # generate the Birkhoff polytope of order 5
#' P = gen_birkhoff(5)
#' @export
gen_birkhoff <- function(n) {
  
  kind_gen = 7
  m_gen = 0
  
  Mat = poly_gen(kind_gen, FALSE, FALSE, n, m_gen)
  
  # first column is the vector b
  b = Mat[,1]
  Mat = Mat[,-c(1)]
  
  P = Hpolytope$new(Mat, b)
  
  return(P)
  
}
