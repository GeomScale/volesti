#' Compute the Chebychev ball of a H-polytope
#' 
#' For a H-polytope described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}, this function computes the largest inscribed ball (Chebychev ball) of that polytope by solving the corresponding linear program.
#'
#' @param A The matrix of the H-polytope.
#' @param b The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets.
#' @return A \eqn{(d+1)}-dimensional vector that containes the Chebychev ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the Chebychev ball.
#' @examples
#' # compute the Chebychev ball of a 2d unit simplex
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' ball_vec = CheBall(A,b)
#' 
#' # compute the Chebychev ball of 10-dimensional cross polytope
#' PolyList = GenCross(10, 'H')
#' ball_vec = CheBall(PolyList$A, PolyList$b)
#' @export
CheBall <- function(A,b){
  
  Mat = -A
  d = dim(Mat)[2] + 1
  m = dim(Mat)[1]
  r = rep(0,d)
  r[1] = m
  r[2] = d
  
  Mat = matrix(cbind(b, Mat), ncol = dim(Mat)[2] + 1)
  Mat = matrix(rbind(r, Mat), ncol = dim(Mat)[2])
  
  ball_only = TRUE
  
  #---------------------#
  verbose = FALSE
  vpoly = FALSE
  Zono = FALSE
  rotate_only = FALSE
  round_only = FALSE
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
  sample_only = FALSE
  numpoints = 0
  variance = 0
  coordinate = TRUE
  rounding = FALSE
  gen_only = FALSE
  Vpoly_gen = FALSE
  kind_gen = -1
  dim_gen = 0
  m_gen = 0
  exact_zono = FALSE
  sam_simplex = FALSE
  sam_can_simplex = FALSE
  sam_arb_simplex = FALSE
  #-------------------#
  
  Mat = vol_R(Mat, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only, 
              rotate_only, ball_only, sample_only, sam_simplex, sam_can_simplex, sam_arb_simplex, 
              numpoints, variance, coordinate, rounding, verbose)
  
  retvec = c(Mat[1,])
  return(retvec)
    
}
