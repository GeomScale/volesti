#' Spectrahedron Class
#'
#' An S4 class representing a spectrahedron defined by a Linear Matrix Inequality (LMI).
#'
#' @slot ptr External pointer to C++ Spectrahedron object
#' @slot dimension Numeric: dimension of the space
#' @slot lmi_matrices List: matrices defining the LMI (if available)
#' @slot interior_point Numeric: coordinates of an interior point (if known)
#'
#' @export
setClass(
  "Spectrahedron",
  slots = list(
    ptr = "externalptr",
    dimension = "numeric",
    lmi_matrices = "list",
    interior_point = "numeric"
  ),
  validity = function(object) {
    if (!is.null(object@ptr) && typeof(object@ptr) != "externalptr") {
      return("ptr must be an external pointer")
    }
    if (length(object@dimension) != 1 || object@dimension <= 0) {
      return("dimension must be positive")
    }
    TRUE
  }
)

#' Create Spectrahedron from SDPA Format File
#'
#' @param sdpa_file Character string: path to SDPA format file
#' @param interior_point Numeric vector (optional): starting interior point
#'
#' @return \code{Spectrahedron} object
#' @export
#'
#' @examples
#' \dontrun{
#'   S <- spectrahedron_from_sdpa("path/to/file.txt")
#' }
spectrahedron_from_sdpa <- function(sdpa_file, interior_point = NULL) {
  
  if (!file.exists(sdpa_file)) {
    stop("File not found: ", sdpa_file)
  }
  
  ptr <- .load_spectrahedron_from_sdpa_cpp(sdpa_file)
  
  if (is.null(ptr) || !inherits(ptr, "externalptr")) {
    stop("Failed to load spectrahedron from file")
  }
  
  dim <- .get_spectrahedron_dimension_cpp(ptr)
  
  new(
    "Spectrahedron",
    ptr = ptr,
    dimension = dim,
    lmi_matrices = list(),
    interior_point = if (!is.null(interior_point)) as.numeric(interior_point) else numeric(0)
  )
}

#' Create Spectrahedron from Matrices
#'
#' @param L0 Numeric matrix: the constant matrix L_0 in the LMI
#' @param L_list List of numeric matrices: matrices L_1, ..., L_n
#'
#' @return \code{Spectrahedron} object
#' @export
#'
#' @examples
#' \dontrun{
#'   L0 <- matrix(c(1, 0, 0, 1), 2, 2)
#'   L1 <- matrix(c(-1, 0, 0, 0), 2, 2)
#'   L2 <- matrix(c(0, 0, 0, -1), 2, 2)
#'   S <- spectrahedron_from_matrices(L0, list(L1, L2))
#' }
spectrahedron_from_matrices <- function(L0, L_list) {
  
  if (!is.matrix(L0) || !is.numeric(L0)) {
    stop("L0 must be a numeric matrix")
  }
  
  if (!is.list(L_list)) {
    stop("L_list must be a list of matrices")
  }
  
  if (nrow(L0) != ncol(L0)) {
    stop("L0 must be square matrix")
  }
  
  n_rows <- nrow(L0)
  n_cols <- ncol(L0)
  
  for (i in seq_along(L_list)) {
    if (!is.matrix(L_list[[i]]) || !is.numeric(L_list[[i]])) {
      stop("L_list[[", i, "]] must be a numeric matrix")
    }
    if (nrow(L_list[[i]]) != n_rows || ncol(L_list[[i]]) != n_cols) {
      stop("All matrices must have same dimensions as L0")
    }
  }
  
  n_vars <- length(L_list)
  
  ptr <- .create_spectrahedron_from_matrices_cpp(L0, L_list)
  
  if (is.null(ptr) || !inherits(ptr, "externalptr")) {
    stop("Failed to create spectrahedron from matrices")
  }
  
  new(
    "Spectrahedron",
    ptr = ptr,
    dimension = n_vars,
    lmi_matrices = c(list(L0), L_list),
    interior_point = numeric(0)
  )
}

#' @export
setMethod("dim", "Spectrahedron", function(x) {
  as.integer(x@dimension)
})

#' @export
setMethod("show", "Spectrahedron", function(object) {
  cat("Spectrahedron object\n")
  cat("  Dimension:", object@dimension, "\n")
  cat("  Interior point:", if (length(object@interior_point) > 0) "computed" else "not set", "\n")
})
