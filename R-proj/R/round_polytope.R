#' Apply rounding to a convex polytope (H-polytope, V-polytope or a zonotope)
#' 
#' Given a convex H or V polytope or a zonotope as input this function computes a rounding based on minimum volume enclosing ellipsoid of a pointset.
#' 
#' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
#' @param WalkType Optional. A string that declares the random walk method: a) 'CDHR' for Coordinate Directions Hit-and-Run, b) 'RDHR' for Random Directions Hit-and-Run or c) 'BW' for Ball Walk. The default walk is  'CDHR'.
#' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor}.
#' @param radius Optional. The radius for the ball walk.
#' 
#' @return A list with 2 elements: (a) a polytope of the same class as the input polytope class and (b) the element "round_value" which is the determinant of the square matrix of the linear transformation that was applied on the polytope that is given as input.
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' P = HPolytope$new(A, b)
#' listHpoly = round_polytope(P)
#' 
#' # rotate a V-polytope (3d cube) using Random Directions HnR with step equal to 50
#' P = GenCube(3, 'V')
#' ListVpoly = round_polytope(P, WalkType = 'RDHR', walk_length = 50)
#' 
#' # rotate a 4-dimensional zonotope defined by the Minkowski sum of 8 segments using ball walk with a radius equal to 0.02
#' Z = GenZonotope(4,8)
#' ListZono = round_polytope(Z, WalkType = 'BW', radius = 0.02)
#' @export
round_polytope <- function(P, WalkType = NULL, walk_length = NULL, radius = NULL){
  
  Mat = rounding(P, WalkType, walk_length, radius)
  
  # get first row which has the info for round_value
  r = Mat[c(1),]
  round_value = r[1]
  # remove first row
  Mat = Mat[-c(1),]
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  if (vpoly) {
    PP = list("P" = Vpolytope$new(A), "round_value" = round_value)
  }else if (Zono) {
    PP = list("P" = Zonotope$new(A), "round_value" = round_value)
  } else {
    PP = list("P" = Hpolytope$new(A, b), "round_value" = round_value)
  }
  return(PP)
}
