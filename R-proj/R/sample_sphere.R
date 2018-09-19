#' Function to sample uniformly from a sphere or a ball
#' 
#' This function can be used for uniform sampling from a \eqn{d}-dimensional sphere of radius \eqn{r} or from the interior of the sphere (ball).
#'
#' @param dimension The dimension of the sphere.
#' @param radius The radius of the sphere.
#' @param N The number of the points to sample.
#' @param interior A boolean flag. It has to be TRUE when sampling is from the interior of the sphere.
#' 
#' @return A \eqn{d\times N} matrix that containes, column-wise, the sampled points.
#' 
#' @examples 
#' # sample 1200 random points from a 10-dimensional hypersphere of radius 100
#' points = sample_sphere(dimension = 10, radius = 100, N = 1200)
#' 
#' # sample 1000 random points from a 20-dimensional ball of radius 200
#' points = sample_sphere(dimension = 20, radius = 200, N = 1000, interior = TRUE)
sample_sphere <- function(dimension, radius, N, interior) {
  
  if (missing(dimension)) {
    print('Wrong inputs.. see the documentation.')
    return(0)
  } else {
    dim_gen = dimension
  }
  if (missing(radius)) {
    print('Wrong inputs.. see the documentation.')
    return(0)
  } else {
    delta = radius
  }
  if (missing(interior)) {
    sam_ball = FALSE
    sam_sphere = TRUE
  } else {
    sam_ball = TRUE
    sam_sphere = FALSE
  }
  n = 100
  if (!missing(N)){
    n = N
  }
  
  # set flag for verbose mode
  verbose=FALSE
  
  #---------------------#
  A = matrix(c(0,0))
  round_only = FALSE
  rotate_only = FALSE
  W = 0
  e = 0
  internalpoint = c(0)
  Gaussian = FALSE
  win_len = 0
  NN = 0
  C = 0
  ratio = 0
  frac = 0
  ballwalk = FALSE
  vpoly = FALSE
  Zono = FALSE
  exact_zono = FALSE
  sample_only = FALSE
  var = 0
  coord = TRUE
  rounding = FALSE
  gen_only = FALSE
  Vpoly_gen = FALSE
  kind_gen = -1
  m_gen = 0
  exact_zono = FALSE
  ball_only = FALSE
  sam_simplex = FALSE
  sam_can_simplex = FALSE
  sam_arb_simplex = FALSE
  construct_copula = FALSE
  h1 = c(0)
  h2 = c(0)
  slices = 0
  sliceSimplex = FALSE
  #-------------------#
  
  Matpoints = vol_R(A, W, e, internalpoint, Gaussian, win_len, NN, C, ratio, frac,
                    ballwalk, delta, vpoly, Zono, exact_zono, gen_only, Vpoly_gen,
                    kind_gen, dim_gen, m_gen, round_only, rotate_only, ball_only,
                    sample_only, sam_simplex, sam_can_simplex, sam_arb_simplex,
                    sam_ball, sam_sphere, n, var, construct_copula, h1, h2, slices,
                    sliceSimplex, coord, rounding, verbose)
  
  return(Matpoints)
  
}