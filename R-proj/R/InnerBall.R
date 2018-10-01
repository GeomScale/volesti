#' Compute an inscribed ball of a convex polytope
#' 
#' For a H-polytope described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}, this function computes the largest inscribed ball (Chebychev ball) by solving the corresponding linear program.
#' For a V-polytope \eqn{d+1} vertices that define a full dimensional simplex picked at random and the largest inscribed ball of the simplex is computed.
#' For a zonotope we compute \eqn{\delta} = \eqn{min\{r:\ re_i\in P,\ \forall i \}. Then the ball centered at the origin with radius \eqn{\delta /\sqrt{d}} is an inscribed ball.
#'
#' @param P A convex polytope. It is an object from class (a) HPolytope or (b) VPolytope or (c) Zonotope.
#' 
#' @return A \eqn{d+1}-dimensional vector that describes the inscribed ball. The first \eqn{d} coordinates corresponds to the center of the ball and the last one to the radius.
#' 
#' @examples
#' # compute the Chebychev ball of a 2d unit simplex
#' P = GenSimplex(2,'H')
#' ball_vec = CheBall(P)
#' 
#' # compute the Chebychev ball of 3-dimensional cube in V-representation
#' P = GenCube(3, 'V')
#' ball_vec = CheBall(P)
#' @export
InnerBall <- function(P){
  
  if (!missing(P)) {
    repr = class(P)[1]
    if (repr == "HPolytope") {
      vpoly = FALSE
      Zono = FALSE
    } else if(repr == "VPolytope") {
      vpoly = TRUE
      Zono = FALSE
    } else if(repr == "Zonotope") {
      vpoly = FALSE
      Zono = TRUE
    } else {
      stop("Not a known polytope representation.")
    }
    Mat = P$get_mat()
    dimension = dim(Mat)[2] - 1
    walk_length = 10 + floor( dimension / 10 )
  } else {
    stop("No polytope is given.")
  }
  
  Vec = RInnerBall(Mat, Zono, vpoly)
  return(Vec)
  
}
