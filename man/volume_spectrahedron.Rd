#' Volume Computation for Convex Bodies
#'
#' @description
#' Computes the volume of a convex body, with extended support for Spectrahedron objects
#' in addition to the standard polytope types.
#'
#' @param body,x A convex body object of class Polytope, Zonotope, Spectrahedron, etc.
#' @param method Character string specifying the volume computation method.
#'   - For Polytope/Zonotope: "CB" (cooling balls), "SOB" (sequence of balls), "MM" (multiphase), 
#'     "BiW" (biased walk)
#'   - For Spectrahedron: "CB" or "SOB" (same algorithms, optimized for LMI constraints)
#' @param E Numeric, target relative error for volume estimation. Must be positive.
#'   Smaller values give more accurate estimates but take longer to compute.
#'   Default is 0.1 (10% error).
#' @param random_walk Character string specifying the random walk type for sampling:
#'   - "RDHR": Randomized Diameter Hit and Run
#'   - "CDHR": Coordinate Directions Hit and Run  
#'   - "HMC": Hamiltonian Monte Carlo
#'   - "BILLIARD": Billiard walk
#'   Default is "RDHR".
#' @param n_walk_per_round Positive integer, number of random walk steps per round.
#'   Only relevant for certain methods. Default is 1.
#' @param win_len Positive integer, window length for moving average. Default is 4.
#' @param uniform_random_point_in_segment Logical, whether to use uniform random point 
#'   selection. Default is FALSE.
#' @param parameters List of additional method-specific parameters (for advanced use).
#'
#' @return
#' Numeric value representing the estimated volume of the convex body.
#'
#' @details
#' The volume computation for Spectrahedra uses the same algorithmic framework as for
#' polytopes but leverages the special structure of LMI constraints for efficiency.
#'
#' The cooling balls method progressively shrinks the spectrahedron via scaling while
#' accumulating volume estimates. The sequence of balls method uses a different approach
#' to build up the volume estimate.
#'
#' The random walk type significantly affects both computation time and accuracy.
#' RDHR is generally robust but may be slow; BILLIARD can be faster for certain geometries.
#' HMC provides better mixing for high-dimensional problems.
#'
#' For very tight error bounds (small E), computation time can grow significantly.
#' A value of E = 0.05 typically provides good accuracy for most applications.
#'
#' @examples
#' \dontrun{
#'   # Load a spectrahedron
#'   spec <- load_spectrahedron("path/to/file.sdpa")
#'   
#'   # Compute volume with default parameters
#'   vol <- volume(spec)
#'   
#'   # Compute with specific method and error tolerance
#'   vol <- volume(spec, method = "CB", E = 0.05)
#'   
#'   # Use different random walk for faster computation
#'   vol <- volume(spec, method = "CB", E = 0.1, random_walk = "BILLIARD")
#' }
#'
#' @references
#' Emiris, I.Z., Fisikopoulos, V. (2014). 
#' Practical volume computation of structured convex bodies via a geometric approach.
#' ACM Trans. Math. Softw., 40(2):16.
#'
#' @seealso
#' \code{\link{sample_points}} for sampling from the same convex bodies
#' \code{\link{Spectrahedron-class}} for Spectrahedron objects
#'
#' @export
#' @docType methods
setMethod("volume",
  signature(object = "Spectrahedron"),
  definition = function(x, method = "CB", E = 0.1, random_walk = "RDHR",
                       n_walk_per_round = 1, win_len = 4,
                       uniform_random_point_in_segment = FALSE,
                       parameters = list()) {
    
    # Input validation
    if (!is(x, "Spectrahedron")) {
      stop("volume method accepts Spectrahedron only")
    }
    
    if (is.null(x@ptr)) {
      stop("Invalid Spectrahedron object (null pointer)")
    }
    
    # Validate method parameter
    method <- match.arg(method, c("CB", "SOB"))
    
    # Validate error bound
    if (!is.numeric(E) || length(E) != 1 || E <= 0) {
      stop("E must be a positive numeric value. Use E should be a value in (0, 1)")
    }
    
    # Validate walk type
    walk_type <- match.arg(random_walk, 
                          c("RDHR", "CDHR", "HMC", "BILLIARD"))
    
    # Validate other parameters
    if (!is.numeric(n_walk_per_round) || n_walk_per_round < 1) {
      stop("n_walk_per_round must be a positive integer")
    }
    
    if (!is.numeric(win_len) || win_len < 1) {
      stop("win_len must be a positive integer")
    }
    
    if (!is.logical(uniform_random_point_in_segment)) {
      stop("uniform_random_point_in_segment must be logical")
    }
    
    # Call C++ wrapper function
    vol <- volume_spectrahedra_rcpp(
      x@ptr,
      method = method,
      E = E,
      walk_type = walk_type,
      n_walk_per_round = as.integer(n_walk_per_round),
      win_len = as.integer(win_len),
      uniform_random_point = uniform_random_point_in_segment
    )
    
    if (is.na(vol) || !is.finite(vol) || vol <= 0) {
      stop("Volume computation failed or returned invalid result")
    }
    
    return(vol)
  }
)

#' @rdname volume
#' @export
setMethod("volume",
  signature(object = "Polytope"),
  definition = function(x, method = "CB", E = 0.1, random_walk = "RDHR",
                       n_walk_per_round = 1, win_len = 4,
                       uniform_random_point_in_segment = FALSE,
                       parameters = list()) {
    # Existing polytope volume computation
    # (This is a placeholder showing how the extension integrates)
    # In the actual volesti package, this would call existing implementations
    stop("Call the existing polytope volume method")
  }
)
