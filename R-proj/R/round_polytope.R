#' Apply rounding to a convex polytope (H-polytope, V-polytope or a zonotope)
#' 
#' Given a convex H or V polytope or a zonotope as input this function computes a rounding based on minimum volume enclosing ellipsoid of a pointset.
#' 
#' @param A Only for H-polytopes. The \eqn{m\times d} matrix \eqn{A} that containes the directions of the \eqn{m} facets.
#' @param b Only for H-polytopes. The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets s.t.: \eqn{Ax\leq b}.
#' @param V Only for V-polytopes. The \eqn{m\times d} matrix V that containes row-wise the \eqn{m} \eqn{d}-dimensional vertices of the polytope.
#' @param G Only for zonotopes. The \eqn{m\times d} matrix G that containes row-wise the \eqn{m} \eqn{d}-dimensional segments that define a zonotope.
#' @param walk_length Optional. The number of the steps for the random walk. Default value is \eqn{\lfloor 10+d/10\rfloor}.
#' @param ball_walk Optional. Boolean parameter to use ball walk, only for CG algorithm. Default value is false.
#' @param delta Optional. The radius for the ball walk.
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is true.
#' 
#' @return For H-polytopes the return value is a list that containes a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b s.t.: \eqn{Ax\leq b}. For V-polytopes and zonotopes the return value is a list with first element a \eqn{m\times d} matrix that containes row-wise the \eqn{d}-dimensional vertices or segments respectively. For all the representations the returned list containes element "round_value" which is the determinant of the square matrix of the linear transformation that was applied on the polytope that is given as input.
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' listHpoly = round_polytope(A=A, b=b)
#' 
#' # rotate a V-polytope (3d cube) using Random Directions HnR
#' Vmat = GenCube(3, 'V')
#' ListVpoly = round_polytope(V=Vmat, coordinate=FALSE)
#' 
#' # rotate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' Zmat = GenZonotope(10,20)
#' ListZono = round_polytope(G=Zmat)
#' @export
round_polytope <- function(P, method) {
  
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
  
  coordinate = TRUE
  ball_walk = FALSE
  delta = -1
  if(!missing(method)) {
    if(!is.null(method$coord)){
      coordinate = method$coord
    }
    if(!is.null(method$WalkT)) {
      if(method$WalkT=="hnr") {
        ball_walk = FALSE
        delta = -1
      } else if(method$WalkT=="bw") {
        if(!is.null(method$coord)){
          warning("Ball walk and coordinate are both declared. Ball walk is going to be used.")
        }
        coordinate = TRUE
        ball_walk = TRUE
        delta = -1
        if(!is.null(method$delta)){
          delta = method$delta
        }
      } else {
        stop("Not a known random walk method.")
      }
    }
    if(!is.null(method$W)) {
      walk_length = method$W
      if (walk_length<=0) {
        stop("Walk length must be a positive value.")
      }
    }
  }
  
  Mat = rounding(Mat, walk_length, coordinate, ball_walk, delta, vpoly, Zono)
  
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
    PP = list("P" = VPolytope(V=A), "round_value" = round_value)
  }else if (Zono) {
    PP = list("P" = Zonotope(G=A), "round_value" = round_value)
  } else {
    PP = list("P" = HPolytope("A"=A, "b"=b), "round_value" = round_value)
  }
  return(PP)
}
