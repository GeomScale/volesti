#' Compute the exact volume of a zonotope
#' 
#' Given the \eqn{m \times d} matrix that containes the \eqn{m} segments that define the \eqn{d}-dimensional zonotope, this function computes the sum of the absolute values of the determinants of all the \eqn{d \times d} submatrices.
#' 
#' @param ZonoMat The \eqn{m \times d} matrix that containes the segments that define the zonotope.
#' 
#' @return The exact volume of the zonotope
#' @examples
#' 
#' # compute the exact volume of a 5-dimensional zonotope defined by the Minkowski sum of 10 segments
#' ZonoMat = GenZonotope(5, 10)
#' vol = ExactZonoVol(ZonoMat)
#' @export
exact_vol <- function(Z, exact) {
  
  exact_zono = FALSE
  exact_cube = FALSE
  exact_simplex = FALSE
  exact_cross = FALSE
  dim = 0
  if(!missing(Z)) {
    if (!is.null(exact$cube, exact$cross)) {
      warning("If a polytope is given then cube and cross in exact list should be NULL.")
    }
    Mat = Z$get_mat()
    dimension = dim(Mat)[2] - 1
    if (missing(exact)) {
      exact_zono = TRUE
    } else {
      if(!is.null(exact$simplex)) {
        if(exact$simplex) {
          if(class(Z)[1]!="Volytope") {
            stop("To sample uniformly points from a simplex, it has to be given in V-representation.")
          } else if(dimension != dim(Mat)[1]-2) {
            stop("The polytope is not a simplex.")
          }
          dim = dimension
          exact_simplex = TRUE
        } else {
          exact_zono = TRUE
        }
      }
    }
  } else if (!missing(exact)) {
    if(is.null(exact$dim)) {
      stop("You have to declare dimension for an exact volume computation of a known polytope.")
    } else {
      dim = exact$dim
      if (dim<2) {
        stop("Dimension has to be greater than 2.")
      }
    }
    Mat = matrix(c(0))
    if (!is.null(exact$simplex)) {
      exact_simplex = exact$simplex
    }
    if (!is.null(exact$cube)) {
      exact_cube = exact$cube
      if (exact_cube && exact_simplex) {
        stop("Wrong declaration.. Only exact volume of one polytope can be computed.")
      }
    }
    if (!is.null(exact$cross)) {
      exact_cross = exact$cross
      if (exact_cross && (exact_cube || exact_simplex)) {
        stop("Wrong declaration.. Only exact volume of one polytope can be computed.")
      }
    }
    if (!exact_cross && !exact_cube && !exact_simplex) {
      stop("Np exact volume computation is enabled.")
    }
  } else {
    stop("Wrong inputs..See the documentation.")
  }
  
  vol = Rvol_exact(Mat, exact_zono, exact_cube, exact_simplex, exact_cross, dim)
  return(vol)
}