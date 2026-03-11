#' Spectrahedron S4 Class
#'
#' An S4 class representing a spectrahedron (spectrahedra in plural), 
#' which is a convex body defined by a Linear Matrix Inequality (LMI) constraint.
#'
#' A spectrahedron is defined by the set: \deqn{\{x \in R^n | L_0 + \sum_{i=1}^{n} x_i L_i \succeq 0\}}{
#' {x in R^n | L_0 + sum(x_i * L_i) >= 0 (positive semidefinite)}}
#' where \eqn{L_0, L_1, ..., L_n} are symmetric matrices and \eqn{\succeq 0}{>= 0} denotes 
#' positive semidefiniteness.
#'
#' @slot ptr An external pointer to the underlying C++ Spectrahedron object
#' @slot dimension The dimension of the space in which the spectrahedron is embedded
#'
#' @details
#' Spectrahedra form a fundamental class of convex bodies that unify linear inequalities
#' (which define polytopes) with polynomial constraints. They are the feasible regions
#' of semidefinite programs (SDPs).
#'
#' Key properties:
#' - Spectrahedra are non-polyhedral in general
#' - They contain polytopes as special cases (when all matrices are diagonal)
#' - They can be sampled uniformly using random walks
#' - Volume can be computed using cooling balls or sequence of balls methods
#'
#' @examples
#' \dontrun{
#'   # Load a spectrahedron from SDPA format
#'   spec <- load_spectrahedron("path/to/file.sdpa")
#'   
#'   # Get dimension
#'   d <- dim(spec)
#'   
#'   # Sample points uniformly
#'   samples <- sample_points(spec, N = 100, walk_type = "RDHR")
#'   
#'   # Compute volume
#'   vol <- volume(spec, method = "CB")
#'   
#'   # Show spectrahedron information
#'   show(spec)
#' }
#'
#' @seealso 
#' \code{\link{load_spectrahedron}} for loading from file
#' \code{\link{spectrahedron_from_matrices}} for direct construction
#' \code{\link{volume,Spectrahedron-method}} for volume computation
#' \code{\link{sample_points,Spectrahedron-method}} for sampling
#'
#' @export
setClass("Spectrahedron",
  slots = list(
    ptr = "externalptr",
    dimension = "integer"
  ),
  validity = function(object) {
    if (!is.null(object@ptr) && is.null(attr(object@ptr, "tag"))) {
      return("Invalid external pointer")
    }
    if (object@dimension <= 0) {
      return("Dimension must be positive")
    }
    return(TRUE)
  }
)

#' Load a Spectrahedron from SDPA Format File
#'
#' @description
#' Creates a Spectrahedron object by loading an SDPA (SemiDefinite Programming format) file.
#'
#' @param filename Character string specifying the path to the SDPA format file.
#'   The file should define a linear matrix inequality constraint.
#'
#' @return
#' A \code{Spectrahedron} object wrapping the loaded spectrahedron.
#'
#' @details
#' SDPA format is a standard format for specifying semidefinite programs. 
#' The format includes:
#' - Number of variables and dimension of matrices
#' - Objective function coefficients
#' - Linear matrix inequality components
#'
#' For detailed SDPA format specification, see the official SDPA documentation.
#'
#' @examples
#' \dontrun{
#'   spec <- load_spectrahedron("path/to/problem.sdpa")
#'   print(spec)
#' }
#'
#' @seealso 
#' \code{\link{spectrahedron_from_matrices}} for creating from R matrices
#'
#' @export
load_spectrahedron <- function(filename) {
  if (!is.character(filename) || length(filename) != 1) {
    stop("filename must be a single character string")
  }
  
  if (!file.exists(filename)) {
    stop("File not found: ", filename)
  }
  
  # Call C++ function to load
  ptr <- load_spectrahedron_from_sdpa_cpp(filename)
  
  # Get dimension
  d <- get_spectrahedron_dimension_cpp(ptr)
  
  if (d <= 0) {
    stop("Failed to load spectrahedron from file")
  }
  
  # Create and return Spectrahedron object
  new("Spectrahedron", ptr = ptr, dimension = as.integer(d))
}

#' Create a Spectrahedron from Matrices
#'
#' @description
#' Creates a Spectrahedron object from explicitly specified LMI matrices.
#'
#' @param L0 Numeric matrix defining the constant term \eqn{L_0} in the LMI.
#' @param L_list List of numeric matrices \eqn{L_1, L_2, ..., L_n} in the LMI.
#'   The list length determines the dimension of the spectrahedron.
#'
#' @return
#' A \code{Spectrahedron} object representing the convex body.
#'
#' @details
#' The spectrahedron is defined as:
#' \deqn{\{x \in R^n | L_0 + \sum_{i=1}^{n} x_i L_i \succeq 0\}}{
#' {x in R^n | L_0 + sum(x_i * L_i) >= 0 (positive semidefinite)}}
#'
#' All matrices must be square and have the same dimensions.
#' The matrices should be symmetric (or nearly symmetric up to numerical precision).
#'
#' @examples
#' \dontrun{
#'   # Create a simple 2D spectrahedron
#'   L0 <- matrix(c(1, 0, 0, 1), nrow = 2)
#'   L1 <- matrix(c(0, 1, 1, 0), nrow = 2)
#'   L2 <- matrix(c(1, 0, 0, -1), nrow = 2)
#'   
#'   spec <- spectrahedron_from_matrices(L0, list(L1, L2))
#'   print(spec)
#' }
#'
#' @seealso 
#' \code{\link{load_spectrahedron}} for loading from SDPA file
#'
#' @export
spectrahedron_from_matrices <- function(L0, L_list) {
  if (!is.matrix(L0)) {
    stop("L0 must be a matrix")
  }
  
  if (!is.list(L_list)) {
    stop("L_list must be a list of matrices")
  }
  
  if (nrow(L0) != ncol(L0)) {
    stop("L0 must be square")
  }
  
  m <- nrow(L0)
  
  for (i in seq_along(L_list)) {
    if (!is.matrix(L_list[[i]])) {
      stop("L_list[[", i, "]] is not a matrix")
    }
    if (nrow(L_list[[i]]) != m || ncol(L_list[[i]]) != m) {
      stop("L_list[[", i, "]] has incompatible dimensions")
    }
  }
  
  n <- length(L_list)
  
  # Call C++ function to create
  ptr <- create_spectrahedron_from_matrices_cpp(L0, L_list)
  
  # Create and return Spectrahedron object
  new("Spectrahedron", ptr = ptr, dimension = as.integer(n))
}

#' Dimension of a Spectrahedron
#'
#' @description
#' Returns the dimension (number of variables) of a spectrahedron.
#'
#' @param x A \code{Spectrahedron} object.
#'
#' @return
#' An integer representing the dimension of the spectrahedron.
#'
#' @examples
#' \dontrun{
#'   spec <- load_spectrahedron("path/to/file.sdpa")
#'   d <- dim(spec)
#'   cat("Spectrahedron dimension:", d, "\n")
#' }
#'
#' @export
setMethod("dim", "Spectrahedron", function(x) {
  if (is.null(x@ptr)) {
    stop("Invalid Spectrahedron object (null pointer)")
  }
  return(x@dimension)
})

#' Display a Spectrahedron
#'
#' @description
#' Prints information about a Spectrahedron object.
#'
#' @param object A \code{Spectrahedron} object.
#'
#' @details
#' Displays the dimension and validity of the spectrahedron.
#'
#' @export
setMethod("show", "Spectrahedron", function(object) {
  cat("Spectrahedron (n =", object@dimension, ")\n")
  if (!is.null(object@ptr)) {
    cat("Valid: Yes\n")
  } else {
    cat("Valid: No (null pointer)\n")
  }
  cat("Methods available: volume, sample_points, dim\n")
})
