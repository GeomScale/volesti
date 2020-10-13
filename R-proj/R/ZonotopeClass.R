#' An R class to represent a Zonotope
#'
#' A zonotope is a convex polytope defined by the Minkowski sum of \eqn{m} \eqn{d}-dimensional segments.
#'
#' \describe{
#'    \item{G}{An \eqn{m\times d} numerical matrix that contains the segments (or generators) row-wise}
#' 
#'    \item{volume}{The volume of the polytope if it is known, \eqn{NaN} otherwise by default.}
#'    
#'    \item{type}{A character with default value 'Zonotope', to declare the representation of the polytope.}
#' }
#'
#' @examples 
#' G = matrix(c(2,3,-1,7,0,0),ncol = 2, nrow = 3, byrow = TRUE)
#' P = Zonotope(G = G)
#'  
#' @name Zonotope-class
#' @rdname Zonotope-class
#' @exportClass Zonotope
Zonotope <- setClass (
  # Class name
  "Zonotope",
  
  # Defining slot type
  representation (
    G = "matrix",
    volume = "numeric",
    type = "character"
  ),
  
  # Initializing slots
  prototype = list(
    G = as.matrix(0),
    volume = as.numeric(NaN),
    type = "Zonotope"
  )
)
