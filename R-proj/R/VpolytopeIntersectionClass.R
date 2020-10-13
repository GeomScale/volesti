#' An R class to represent the intersection of two V-polytopes
#'
#' An intersection of two V-polytopes is defined by the intersection of the two coresponding convex hulls.
#'
#' \describe{
#'    \item{V1}{An \eqn{m\times d} numerical matrix that contains the vertices of the first V-polytope (row-wise).}
#'    
#'    \item{V2}{An \eqn{q\times d} numerical matrix that contains the vertices of the second V-polytope (row-wise).}
#' 
#'    \item{volume}{The volume of the polytope if it is known, \eqn{NaN} otherwise by default.}
#'    
#'    \item{type}{A character with default value 'VpolytopeIntersection', to declare the representation of the polytope.}
#' }
#' 
#' @examples
#' P1 = gen_simplex(2,'V')
#' P2 = gen_cross(2,'V')
#' P = VpolytopeIntersection(V1 = P1@V, V2 = P2@V)
#' 
#' @name VpolytopeIntersection-class
#' @rdname VpolytopeIntersection-class
#' @exportClass VpolytopeIntersection
VpolytopeIntersection <- setClass (
  # Class name
  "VpolytopeIntersection",
  
  # Defining slot type
  representation (
    V1 = "matrix",
    V2 = "matrix",
    volume = "numeric",
    type = "character"
  ),
  
  # Initializing slots
  prototype = list(
    V1 = as.matrix(0),
    V2 = as.matrix(0),
    volume = as.numeric(NaN),
    type = "VpolytopeIntersection"
  )
)
