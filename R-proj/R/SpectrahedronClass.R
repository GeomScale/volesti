#' An R class to represent a Spectrahedron
#'
#' A spectrahedron is a convex body defined by a linear matrix inequality of the form \eqn{A_0 + x_1 A_1 + ... + x_n A_n \preceq 0}.
#' The matrices \eqn{A_i} are symmetric \eqn{m \times m} real matrices and \eqn{\preceq 0} denoted negative semidefiniteness.
#'
#' \describe{
#'    \item{matrices}{A List that contains the matrices \eqn{A_0, A_1, ..., A_n}.}
#' }
#' 
#' @examples
#' A0 = matrix(c(-1,0,0,0,-2,1,0,1,-2), nrow=3, ncol=3, byrow = TRUE)
#' A1 = matrix(c(-1,0,0,0,0,1,0,1,0), nrow=3, ncol=3, byrow = TRUE)
#' A2 = matrix(c(0,0,-1,0,0,0,-1,0,0), nrow=3, ncol=3, byrow = TRUE)
#' lmi = list(A0, A1, A2)
#' S = Spectrahedron(matrices = lmi);
#'  
#' @name Spectrahedron-class
#' @rdname Spectrahedron-class
#' @exportClass Spectrahedron
Spectrahedron <- setClass (
  # Class name
  "Spectrahedron",
  
  # Defining slot type
  representation (
    matrices = "list"
  ),
  
  # Initializing slots
  prototype = list(
    matrices = list()
  )
)
