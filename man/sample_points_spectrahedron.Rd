#' Sampling from Convex Bodies
#'
#' @description
#' Generates uniformly distributed random samples from the interior of a convex body,
#' with extended support for Spectrahedron objects in addition to standard polytope types.
#'
#' @param x A convex body object of class Polytope, Zonotope, Spectrahedron, etc.
#' @param N Positive integer, number of samples to generate.
#' @param walk_type Character string specifying the random walk type:
#'   - "RDHR": Randomized Diameter Hit and Run (default)
#'   - "CDHR": Coordinate Directions Hit and Run
#'   - "HMC": Hamiltonian Monte Carlo
#'   - "BILLIARD": Billiard walk
#' @param n_burns Non-negative integer, number of burn-in steps before collecting samples.
#'   These initial steps are discarded to reduce dependence on the starting point.
#'   Default is 0.
#' @param starting_point Numeric vector, optional starting point for the random walk.
#'   If not provided, an interior point is computed automatically.
#' @param parameters List of additional walk-specific parameters (for advanced use).
#'
#' @return
#' A numeric matrix of dimension D × N where D is the dimension of the convex body
#' and N is the number of samples. Each column is a sample point uniformly drawn
#' from the interior of the convex body.
#'
#' @details
#' The function uses Markov chain Monte Carlo (MCMC) methods to generate approximately
#' uniform samples from convex bodies. Different random walks have different properties:
#'
#' \describe{
#'   \item{RDHR}{Randomized Diameter Hit-and-Run: Choose random direction, walk until
#'     boundary, then pick random point. Generally stable and well-mixing.}
#'   \item{CDHR}{Coordinate Directions Hit-and-Run: Similar to RDHR but uses coordinate
#'     axes instead of random directions. Faster in some cases.}
#'   \item{HMC}{Hamiltonian Monte Carlo: Uses momentum-based moves for better mixing.
#'     Advantage in high dimensions.}
#'   \item{BILLIARD}{Billiard walk: Simulates ball bouncing off boundary. Very fast
#'     but may have lower mixing rate for some geometries.}
#' }
#'
#' For Spectrahedra, all four walk types are supported and leverage the structure
#' of LMI constraints for efficient boundary computations.
#'
#' The burn-in period (n_burns) helps reduce the dependence of samples on the starting
#' point and allows the chain to reach its stationary distribution. Using n_burns = 100
#' or higher is recommended for high-accuracy applications.
#'
#' @examples
#' \dontrun{
#'   # Load a spectrahedron
#'   spec <- load_spectrahedron("path/to/file.sdpa")
#'   
#'   # Generate 100 uniform samples using default parameters
#'   samples <- sample_points(spec, N = 100)
#'   
#'   # Use different walk with burn-in steps
#'   samples <- sample_points(spec, N = 100, walk_type = "HMC", n_burns = 100)
#'   
#'   # Specify a custom starting point
#'   starting_point <- rep(0, dim(spec))
#'   samples <- sample_points(spec, N = 100, starting_point = starting_point)
#' }
#'
#' @references
#' Vempala, S.S. (2005). 
#' Geometric Random Walks: A Survey. MSRI Publications 52.
#'
#' Simonovits, M. (2003).
#' How to compute the volume in high dimension? Mathematical Aspects of Computer Science.
#'
#' @seealso
#' \code{\link{volume}} for computing the volume of the same convex bodies
#' \code{\link{Spectrahedron-class}} for Spectrahedron objects
#'
#' @export
#' @docType methods
setMethod("sample_points",
  signature(object = "Spectrahedron"),
  definition = function(x, N = 100, walk_type = "RDHR", n_burns = 0,
                       starting_point = NULL, parameters = list()) {
    
    # Input validation
    if (!is(x, "Spectrahedron")) {
      stop("sample_points method accepts Spectrahedron only")
    }
    
    if (is.null(x@ptr)) {
      stop("Invalid Spectrahedron object (null pointer)")
    }
    
    # Validate N parameter
    if (!is.numeric(N) || length(N) != 1 || N < 1 || N != floor(N)) {
      stop("N must be a positive integer")
    }
    N <- as.integer(N)
    
    # Validate walk_type
    walk_type <- match.arg(walk_type, 
                          c("RDHR", "CDHR", "HMC", "BILLIARD"))
    
    # Validate burn-in
    if (!is.numeric(n_burns) || length(n_burns) != 1 || n_burns < 0 || n_burns != floor(n_burns)) {
      stop("n_burns must be a non-negative integer")
    }
    n_burns <- as.integer(n_burns)
    
    # Validate starting point if provided
    if (!is.null(starting_point)) {
      starting_point <- as.numeric(starting_point)
      if (length(starting_point) != x@dimension) {
        stop("starting_point dimension (", length(starting_point),") does not match ",
             "spectrahedron dimension (", x@dimension, ")")
      }
    }
    
    # Call C++ wrapper function
    samples_list <- sample_spectrahedra_rcpp(
      x@ptr,
      N = N,
      walk_type = walk_type,
      n_burns = n_burns,
      starting_point = starting_point
    )
    
    if (is.null(samples_list) || !is.list(samples_list)) {
      stop("Sampling failed")
    }
    
    # Convert list to matrix format
    # Expected structure: list with samples and related info
    if (!is.null(samples_list$samples)) {
      samples_matrix <- samples_list$samples
    } else {
      samples_matrix <- as.matrix(samples_list)
    }
    
    # Verify output dimensions
    if (nrow(samples_matrix) != x@dimension || ncol(samples_matrix) != N) {
      stop("Output dimensions do not match expected (", x@dimension, " x ", N, ")")
    }
    
    # Verify samples are finite
    if (any(!is.finite(samples_matrix))) {
      stop("Sampling produced non-finite values")
    }
    
    return(samples_matrix)
  }
)

#' @rdname sample_points 
#' @export
setMethod("sample_points",
  signature(object = "Polytope"),
  definition = function(x, N = 100, walk_type = "RDHR", n_burns = 0,
                       starting_point = NULL, parameters = list()) {
    # Existing polytope sampling
    # (This is a placeholder showing how the extension integrates)
    stop("Call the existing polytope sampling method")
  }
)
