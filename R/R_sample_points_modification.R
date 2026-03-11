#' @rdname sample_points
#' @description
#'   Sample random points from a Spectrahedron using various random walk methods
#'
#' @param object A \code{Spectrahedron} object
#' @param n Integer: number of points to sample
#' @param walk_type Character: type of random walk ("RDHR", "HMC", "CDHR", "BILLIARD")
#' @param walk_length Integer: length of each walk segment
#' @param random_seed Integer: seed for random number generator
#' @param starting_point Numeric vector: initial point (auto-computed if NULL)
#' @param n_burns Integer: number of burn-in samples
#' @param ... Additional arguments (ignored)
#'
#' @return Matrix: rows are sampled points (n x d matrix where d is dimension)
#'
#' @export
setMethod(
  "sample_points",
  signature(object = "Spectrahedron"),
  function(object,
           n,
           walk_type = "RDHR",
           walk_length = 10L,
           random_seed = 1L,
           starting_point = NULL,
           n_burns = 0L,
           ...) {
    
    if (object@dimension == 0 || is.null(object@ptr)) {
      stop("Invalid Spectrahedron object")
    }
    
    n <- as.integer(n)
    if (n <= 0L) {
      stop("n must be a positive integer")
    }
    
    walk_type <- match.arg(
      walk_type,
      c("RDHR", "HMC", "CDHR", "BILLIARD"),
      several.ok = FALSE
    )
    
    walk_length <- as.integer(walk_length)
    if (walk_length < 1L) {
      stop("walk_length must be >= 1")
    }
    
    random_seed <- as.integer(random_seed)
    n_burns <- as.integer(n_burns)
    if (n_burns < 0L) {
      stop("n_burns must be >= 0")
    }
    
    if (!is.null(starting_point)) {
      starting_point <- as.numeric(starting_point)
      if (length(starting_point) != object@dimension) {
        stop("starting_point dimension must match Spectrahedron dimension")
      }
    }
    
    tryCatch({
      sample_matrix <- sample_spectrahedra_rcpp(
        object@ptr,
        n = n,
        walk_type = walk_type,
        walk_length = walk_length,
        seed = random_seed,
        starting_point = starting_point,
        n_burns = n_burns
      )
      
      if (!is.matrix(sample_matrix)) {
        stop("C++ function returned invalid type")
      }
      
      if (nrow(sample_matrix) != object@dimension) {
        sample_matrix <- t(sample_matrix)
      }
      
      return(sample_matrix)
      
    }, error = function(e) {
      stop("Error sampling points: ", conditionMessage(e))
    })
  }
)
