#' Compute the percentage of the volume of the unit simplex that is contained in the intersection of a half-space and the unit simplex
#'
#' When a half-space \eqn{H} is given as a pair of a vector \eqn{c\in R^d} and a scalar \eqn{z0\in R} s.t.: \eqn{c^Tx\leq z0} this function calls the Ali's version of the Varsi formula.
#' 
#' @param H A \eqn{d}-dimensional vector that defines the direction of the hyperplane.
#' @param z0 The scalar that defines the half-space.
#' 
#' @references \cite{Varsi, Giulio,
#' \dQuote{The multidimensional content of the frustum of the simplex,} \emph{Pacific J. Math. 46, no. 1, 303--314,} 1973.}
#' 
#' @references \cite{Ali, Mir M.,
#' \dQuote{Content of the frustum of a simplex,} \emph{ Pacific J. Math. 48, no. 2, 313--322,} 1973.}
#' 
#' @return The percentage of the volume of the unit simplex that is contained in the intersection of the given half-space and the unit simplex
#' 
#' @examples 
#' # compute the frustum of H: -x+y<=0
#' H=c(-1,1)
#' z0=0
#' frustum = sliceOfSimplex(H=H, z0=z0)
sliceOfSimplex <- function(H, z0) {
  
  H = c(H, z0)
  
  sliceSimplex = TRUE
  
  # set flag for verbose mode
  verbose=FALSE
  
  #---------------------#
  construct_copula = FALSE
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
  dim_gen = 0
  exact_zono = FALSE
  ball_only = FALSE
  sam_simplex = FALSE
  sam_can_simplex = FALSE
  sam_arb_simplex = FALSE
  sam_ball = FALSE
  sam_sphere = FALSE
  h2 = c(0)
  E = matrix(c(0,0))
  slices = 0
  n = 0
  #-------------------#
  
  vol = vol_R(E, W, e, internalpoint, Gaussian, win_len, NN, C, ratio, frac,
              ballwalk, delta, vpoly, Zono, exact_zono, gen_only, Vpoly_gen,
              kind_gen, dim_gen, m_gen, round_only, rotate_only, ball_only,
              sample_only, sam_simplex, sam_can_simplex, sam_arb_simplex, sam_ball,
              sam_sphere, n, var, construct_copula, H, h2, slices, sliceSimplex,
              coord, rounding, verbose)
  
  return(vol[1,1])
}