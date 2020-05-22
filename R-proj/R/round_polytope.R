#' Apply rounding to a convex polytope (H-polytope, V-polytope or a zonotope)
#' 
#' Given a convex H or V polytope or a zonotope as input this functionbrings the polytope in well rounded position based on minimum volume enclosing ellipsoid of a pointset.
#' 
#' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
#' @param seed Optional. A fixed seed for the number generator.
#' 
#' @return A list with 2 elements: (a) a polytope of the same class as the input polytope class and (b) the element "round_value" which is the determinant of the square matrix of the linear transformation that was applied on the polytope that is given as input.
#'
#' @references \cite{Michael J. Todd and E. Alper Yildirim,
#' \dQuote{On Khachiyan’s Algorithm for the Computation of Minimum Volume Enclosing Ellipsoids,} \emph{Discrete Applied Mathematics,} 2007.}
#'
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' P = Hpolytope$new(A, b)
#' listHpoly = round_polytope(P)
#' 
#' # rotate a V-polytope (3d unit cube) using Random Directions HnR with step equal to 50
#' P = gen_cube(3, 'V')
#' ListVpoly = round_polytope(P)
#' 
#' # round a 2-dimensional zonotope defined by 6 generators using ball walk
#' Z = gen_rand_zonotope(2,6)
#' ListZono = round_polytope(Z)
#' @export
round_polytope <- function(P, seed = NULL){
  
  ret_list = rounding(P)
  
  #get the matrix that describes the polytope
  Mat = ret_list$Mat
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  
  type = P$type
  if (type == 2) {
    PP = list("P" = Vpolytope$new(A), "T" = ret_list$T, "shift" = ret_list$shift, "round_value" = ret_list$round_value)
  }else if (type == 3) {
    PP = list("P" = Zonotope$new(A), "T" = ret_list$T, "shift" = ret_list$shift, "round_value" = ret_list$round_value)
  } else {
    PP = list("P" = Hpolytope$new(A,b), "T" = ret_list$T, "shift" = ret_list$shift, "round_value" = ret_list$round_value)
  }
  return(PP)
}
