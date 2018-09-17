#' Function to sample from the unit or an arbitrary simplex
#' 
#' The \eqn{d}-dimensional unit simplex is the set of points \eqn{\vec{x}\in \R}, s.t.: \eqn{\sum_i x_i\leq 1}. The \eqn{d}-dimensional canonical simplex is the set of points \eqn{\vec{x}\in \R}, s.t.: \eqn{\sum_i x_i = 1}.
#' 
#' @param vertices Only for an arbitrary simplex. A \eqn{(d+1)\times d} matrix that containes the vertices of a \eqn{d}-dimensional simplex.
#' @param dimension The dimension of the unit or the canonical simplex.
#' @param N The number of points to sample. Default value is \eqn{100}.
#' @param canonical A boolean flag. It has to be TRUE when sampling is from the \eqn{d}-dimensional canonical simplex. Default value is FALSE.
#' 
#' @description This function can be used to sample uniform points from the \eqn{d}-dimensional unit simplex \eqn{S\subset\R^d} or the \eqn{d-1}-dimensional canonical simplex \eqn{S\subset\R^{d}}.
#' Moreover it can be used to sample uniform points from an arbitrary simplex when \eqn{d+1} \eqn{d}-dimensional vertices that define a full dimensional simplex are given.
#' 
#' @return A \eqn{d\times N} matrix that containes, column-wise, the sampled points.
#' 
#' @examples 
#' # sample 1000 points from the 10-dimensional unit simplex
#' MatPoints = sample_simplex(dimension = 10, N = 1000)
#' 
#' # sample 2000 points from the 20-dimensional canonical simplex
#' MatPoints = sample_simplex(dimension = 10, N = 2000, canonical = TRUE)
#' 
#' # sample 3000 points from an arbitrary simplex
#' V = matrix(c(0,0,0,7,11,0), ncol=2, nrow=3, byrow=TRUE)
#' MatPoints = sample_simplex(vertices = V, N = 3000)
sample_simplex <- function(vertices, dimension, N, canonical) {
  
  if (!missing(vertices)) {
    sam_simplex = FALSE
    sam_can_simplex = FALSE
    sam_arb_simplex = TRUE
    dim_gen = 0
  } else if(!missing(canonical)) {
    if (missing(dimension)) {
      print('Wrong inputs.. see the documentation.')
      return(0)
    }
    sam_simplex = FALSE
    sam_can_simplex = TRUE
    sam_arb_simplex = FALSE
    vertices = matrix(c(0,0))
    dim_gen = dimension
  } else if(missing(dimension)){
    print('Wrong inputs.. see the documentation.')
    return(0)
  } else {
    vertices = matrix(c(0,0))
    sam_simplex = TRUE
    sam_can_simplex = FALSE
    sam_arb_simplex = FALSE
    dim_gen = dimension
  }
  
  n = 100
  if (!missing(N)){
    n = N
  }
  
  # set flag for verbose mode
  verbose=FALSE
  
  #---------------------#
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
  delta =0
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
  sam_ball = FALSE
  sam_sphere = FALSE
  #-------------------#
  
  
  # set timer
  tim = proc.time()
  
  Matpoints = vol_R(vertices, W, e, internalpoint, Gaussian, win_len, NN, C, ratio, frac,
                 ballwalk, delta, vpoly, Zono, exact_zono, gen_only, Vpoly_gen,
                 kind_gen, dim_gen, m_gen, round_only, rotate_only, ball_only,
                 sample_only, sam_simplex, sam_can_simplex, sam_arb_simplex, sam_ball,
                 sam_sphere, n, var, coord, rounding, verbose)
  
  tim = proc.time() - tim
  if (verbose) {
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  return(Matpoints)
}