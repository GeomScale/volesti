#' An R class to represent a V-polytope
#'
#' A V-polytope is a convex polytope defined by the set of its vertices.
#'
#' \describe{
#'    \item{V}{An \eqn{m\times d} numerical matrix that contains the vertices row-wise.}
#' 
#'    \item{volume}{The volume of the polytope if it is known, \eqn{NaN} otherwise by default.}
#'    
#'    \item{type}{A character with default value 'Vpolytope', to declare the representation of the polytope.}
#' }
#'  
#' @examples
#' V = matrix(c(2,3,-1,7,0,0),ncol = 2, nrow = 3, byrow = TRUE)
#' P = Vpolytope(V = V)
#' 
#' @name Vpolytope-class
#' @rdname Vpolytope-class
#' @exportClass Vpolytope
Vpolytope <- setClass (
  # Class name
  "Vpolytope",
  
  # Defining slot type
  representation (
    V = "matrix",
    volume = "numeric",
    type = "character"
  ),
  
  # Initializing slots
  prototype = list(
    V = as.matrix(0),
    volume = as.numeric(NaN),
    type = "Vpolytope"
  )
)
