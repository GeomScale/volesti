#' An R class to represent an H-polytope
#'
#' An H-polytope is a convex polytope defined by a set of linear inequalities or equivalently a \eqn{d}-dimensional H-polytope with \eqn{m} facets is defined by a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b, s.t.: \eqn{Ax\leq b}.
#'
#' \describe{
#'    \item{A}{An \eqn{m\times d} numerical matrix.}
#'
#'    \item{b}{An \eqn{m}-dimensional vector b.}
#' 
#'    \item{volume}{The volume of the polytope if it is known, \eqn{NaN} otherwise by default.}
#'    
#'    \item{type}{A character with default value 'Hpolytope', to declare the representation of the polytope.}
#' }
#' 
#' @examples
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' P = Hpolytope(A = A, b = b)
#'  
#' @name Hpolytope-class
#' @rdname Hpolytope-class
#' @exportClass Hpolytope
Hpolytope <- setClass (
  # Class name
  "Hpolytope",
  
  # Defining slot type
  representation (
    A = "matrix",
    b = "numeric",
    volume = "numeric",
    type = "character"
    ),
  
  # Initializing slots
  prototype = list(
    A = as.matrix(0),
    b = as.numeric(NULL),
    volume = as.numeric(NaN),
    type = "Hpolytope"
  )
)
