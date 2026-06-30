#' @rdname volume
#' @description
#'   Compute volume of a Spectrahedron using Cooling methods
#'
#' @param object A \code{Spectrahedron} object
#' @param volume_method Character: volume computation algorithm ("CB", "SOB")
#' @param walk_length Integer: length of each random walk segment
#' @param tolerance Numeric: relative error tolerance for volume estimate
#' @param n_threads Integer: number of parallel threads (reserved for future)
#' @param parameters List: additional parameters, e.g. list(seed=123)
#' @param ... Additional arguments (ignored)
#'
#' @return Numeric: estimated volume of the spectrahedron
#'
#' @export
setMethod(
  "volume",
  signature(object = "Spectrahedron"),
  function(object, 
           volume_method = "CB",
           walk_length = 10L,
           tolerance = 0.1,
           n_threads = 1L,
           parameters = list(),
           ...) {
    
    if (object@dimension == 0 || is.null(object@ptr)) {
      stop("Invalid Spectrahedron object")
    }
    
    volume_method <- match.arg(
      volume_method,
      c("CB", "SOB"),
      several.ok = FALSE
    )
    
    walk_length <- as.integer(walk_length)
    if (walk_length < 1L) {
      stop("walk_length must be >= 1")
    }
    
    if (!is.numeric(tolerance) || tolerance <= 0) {
      stop("tolerance must be a positive number")
    }
    
    seed <- as.integer(parameters$seed %||% 1L)
    
    tryCatch({
      vol <- volume_spectrahedra_rcpp(
        object@ptr,
        walk_length = walk_length,
        tolerance = tolerance,
        seed = seed,
        volume_method = volume_method
      )
      
      if (!is.numeric(vol) || !is.finite(vol)) {
        warning("Volume computation resulted in invalid value")
        return(NA_real_)
      }
      
      return(vol)
      
    }, error = function(e) {
      stop("Error computing volume: ", conditionMessage(e))
    })
  }
)

# Helper function
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
