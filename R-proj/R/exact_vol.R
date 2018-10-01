#' Compute the exact volume of (a) a zonotope (b) an arbitrary simplex (c) a unit simplex (d) a cross polytope (e) a hypercube
#' 
#' Given a zonotope (as an object of class Zonotope), this function computes the sum of the absolute values of the determinants of all the \eqn{d \times d} submatrices of the \eqn{m\times d} matrix \eqn{G} that contains row-wise the segments that define the zonotope.
#' For an arbitrary simplex that is given in V-representation this function computes the absolute value of the determinant formed by the simplex's points assuming it is shifted to the origin.
#' For a \eqn{d}-dimensional unit simplex, hypercube or cross polytope this function computes the exact well known formulas.
#' 
#' @param Z An object of class Zonotope.
#' @param exact A list that contains parameters for the exact volume computations. When a zonotope is given it should be null.
#' \itemize{
#'  \item{simplex }{A boolean parameter. It has to be TRUE when a simplex is given in V-representation or in order to compute the exact volume of a unit simplex.}
#' }
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
    if (!missing(exact)) {
      if (!is.null(exact$cube, exact$cross)) {
        warning("If a polytope is given then cube and cross in exact list should be NULL.")
      }
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