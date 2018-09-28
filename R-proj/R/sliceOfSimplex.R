#' Compute the percentage of the volume of the unit simplex that is contained in the intersection of a half-space and the unit simplex
#'
#' When a half-space \eqn{H} is given as a pair of a vector \eqn{c\in R^d} and a scalar \eqn{z0\in R} s.t.: \eqn{c^Tx\leq z0} this function calls the Ali's version of the Varsi formula.
#' 
#' @param c A \eqn{d}-dimensional vector that defines the direction of the hyperplane.
#' @param z0 The scalar that defines the half-space.
#' 
#' @references \cite{Varsi, Giulio,
#' \dQuote{The multidimensional content of the frustum of the simplex,} \emph{Pacific J. Math. 46, no. 1, 303--314,} 1973.}
#' 
#' @references \cite{Ali, Mir M.,
#' \dQuote{Content of the frustum of a simplex,} \emph{ Pacific J. Math. 48, no. 2, 313--322,} 1973.}
#' 
#' @return The percentage of the volume of the unit simplex that is contained in the intersection of the given half-space and the unit simplex
#' 
#' @examples 
#' # compute the frustum of H: -x+y<=0
#' c=c(-1,1)
#' z0=0
#' frustum = sliceOfSimplex(c, z0)
sliceOfSimplex <- function(c, z0) {
  
  H = c(c, z0)
  
  vol = SliceSimplex(H)
  
  return(vol)
}
