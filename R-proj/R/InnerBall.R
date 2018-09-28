#' Compute the Chebychev ball of a H-polytope
#' 
#' For a H-polytope described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}, this function computes the largest inscribed ball (Chebychev ball) of that polytope by solving the corresponding linear program.
#'
#' @param A The matrix of the H-polytope.
#' @param b The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets.
#' @return A \eqn{(d+1)}-dimensional vector that containes the Chebychev ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the Chebychev ball.
#' @examples
#' # compute the Chebychev ball of a 2d unit simplex
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' ball_vec = CheBall(A,b)
#' 
#' # compute the Chebychev ball of 10-dimensional cross polytope
#' PolyList = GenCross(10, 'H')
#' ball_vec = CheBall(PolyList$A, PolyList$b)
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
