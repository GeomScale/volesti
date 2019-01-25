#' Compute an inscribed ball of a convex polytope
#' 
#' For a H-polytope described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}, this function computes the largest inscribed ball (Chebychev ball) by solving the corresponding linear program.
#' For a V-polytope \eqn{d+1} vertices that define a full dimensional simplex picked at random and the largest inscribed ball of the simplex is computed.
#' For a zonotope \eqn{P} we compute the minimum \eqn{r} s.t.: \eqn{ r e_i \in P} for all \eqn{i=1, \dots ,d}. Then the ball centered at the origin with radius \eqn{r/ \sqrt{d}} is an inscribed ball.
#'
#' @param P A convex polytope. It is an object from class (a) HPolytope or (b) VPolytope or (c) Zonotope.
#' 
#' @return A \eqn{d+1}-dimensional vector that describes the inscribed ball. The first \eqn{d} coordinates corresponds to the center of the ball and the last one to the radius.
#' 
#' @examples
#' # compute the Chebychev ball of a 2d unit simplex
#' P = GenSimplex(2,'H')
#' ball_vec = InnerBall(P)
#' 
#' # compute the Chebychev ball of 3-dimensional cube in V-representation
#' P = GenCube(3, 'V')
#' ball_vec = InnerBall(P)
#' @export
RRInnerBall <- function(P){
  
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
