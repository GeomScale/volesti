#' Compute the exact volume of a zonotope
#' 
#' Given the \eqn{m \times d} matrix that containes the \eqn{m} segments that define the \eqn{d}-dimensional zonotope, this function computes the sum of the absolute values of the determinants of all the \eqn{d \times d} submatrices.
#' 
#' @param ZonoMat The \eqn{m \times d} matrix that containes the segments that define the zonotope.
#' 
#' @return The exact volume of the zonotope
#' @examples
#' 
#' # compute the exact volume of a 5-dimensional zonotope defined by the Minkowski sum of 10 segments
#' ZonoMat = GenZonotope(5, 10)
#' vol = ExactZonoVol(ZonoMat)
#' @export
ExactZonoVol <- function(ZonoMat) {
  A = ZonoMat
  d = dim(A)[2] + 1
  m = dim(A)[1]
  b = rep(1, m)
  r = rep(0, d)
  r[1] = m
  r[2] = d
  
  A = matrix(cbind(b,A), ncol = dim(A)[2] + 1)
  A = matrix(rbind(r,A), ncol = dim(A)[2])
  
  Zono = TRUE
  exact_zono = TRUE
  
  #---------------------#
  gen_only = FALSE
  round_only = FALSE
  rotate_only = FALSE
  W = 0
  e = 0
  Cheb_ball = c(0)
  annealing = FALSE
  win_len = 0
  N = 0
  C = 0
  ratio = 0
  frac = 0
  ball_walk = FALSE
  delta =0
  Vpoly = FALSE
  Vpoly_gen = FALSE
  kind_gen = -1
  dim_gen = 0
  m_gen = 0
  sample_only = FALSE
  numpoints = 0
  variance = 0
  coordinate = TRUE
  rounding = FALSE
  verbose=FALSE
  ball_only = FALSE
  #-------------------#
  
  vol = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              Vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only, 
              rotate_only, ball_only, sample_only, numpoints, variance, coordinate, rounding, verbose)
  return(vol[1,1])
}