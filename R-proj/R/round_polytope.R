#' Apply rounding to a convex polytope (H-polytope, V-polytope or a zonotope)
#' 
#' Given a convex H or V polytope or a zonotope as input this functionbrings the polytope in well rounded position based on minimum volume enclosing ellipsoid of a pointset.
#' 
#' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
#' @param random_walk Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run, c) \code{'BaW'} for Ball Walk or d) \code{'BiW'} for Billiard walk. The default walk is \code{'CDHR'} for H-polytope and \code{'BiW'} for the other representations.
#' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{1} for \code{'BiW'} and \eqn{\lfloor 10 + d/10\rfloor} otherwise.
#' @param parameters Optional. A list for the parameters of the methods:
#' \itemize{
#'   \item{\code{BW_rad} }{ The radius for the ball walk.}
#'   \item{\code{diameter} }{ The diameter of the polytope. It is used to set the maximum length of the trajectory in billiard walk.}
#' }
#' 
#' @return A list with 2 elements: (a) a polytope of the same class as the input polytope class and (b) the element "round_value" which is the determinant of the square matrix of the linear transformation that was applied on the polytope that is given as input.
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
#' ListVpoly = round_polytope(P, random_walk = 'RDHR', walk_length = 50)
#' 
#' # round a 2-dimensional zonotope defined by 6 generators using ball walk
#' Z = gen_rand_zonotope(2,6)
#' ListZono = round_polytope(Z, random_walk = 'BW')
#' @export
round_polytope <- function(P, random_walk = NULL, walk_length = NULL, parameters = NULL){
  
  ret_list = rounding(P, random_walk, walk_length, parameters)
  
  #get the matrix that describes the polytope
  Mat = ret_list$Mat
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  type = P$type
  if (type == 2) {
    PP = list("P" = Vpolytope$new(A), "round_value" = ret_list$round_value)
  }else if (type == 3) {
    PP = list("P" = Zonotope$new(A), "round_value" = ret_list$round_value)
  } else if(type == 1) {
    PP = list("P" = Hpolytope$new(A,b), "round_value" = ret_list$round_value)
  } else {
    PP = list("P" = VPolyintersectVPoly$new(V1 = t(Mat[,dim(P$V1)[1]]), V2 = t(Mat[,dim(P$V2)[1]])), "round_value" = ret_list$round_value)
  }
  return(PP)
}
