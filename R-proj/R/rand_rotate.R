#' Apply a random rotation to a convex polytope (H-polytope, V-polytope or a zonotope)
#' 
#' Given a convex H or V polytope or a zonotope as input this function applies a random rotation.
#' 
#' @param A Only for H-polytopes. The \eqn{m\times d} matrix \eqn{A} that containes the directions of the \eqn{m} facets.
#' @param b Only for H-polytopes. The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets s.t.: \eqn{Ax\leq b}.
#' @param V Only for V-polytopes. The \eqn{m\times d} matrix V that containes row-wise the \eqn{m} \eqn{d}-dimensional vertices of the polytope.
#' @param G Only for zonotopes. The \eqn{m\times d} matrix G that containes row-wise the \eqn{m} \eqn{d}-dimensional segments that define a zonotope.
#' 
#' @return A random rotation of the polytope that is given as an input. For H-polytopes the return value is a list that containes a \eqn{m\times d} matrix A and a \eqn{m}-dimensional vector b s.t.: \eqn{Ax\leq b}. For V-polytopes and zonotopes the return value is a \eqn{m\times d} matrix that containes row-wise the \eqn{d}-dimensional vertices or segments respectively.
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' listHpoly = rand_rotate(A=A, b=b)
#' 
#' # rotate a V-polytope (3d cube)
#' Vmat = GenCube(3, 'V')
#' matVpoly = rand_rotate(V=Vmat)
#' 
#' # rotate a 5-dimensional zonotope defined by the Minkowski sum of 15 segments
#' Zmat = GenZonotope(5,15)
#' MatZono = rand_rotate(G=Zmat)
#' @export
rand_rotate <- function(A, b, V, G){
  
  vpoly = FALSE
  Zono = FALSE
  if(missing(b)) {
    if(!missing(V)) {
      Mat = V
      vpoly = TRUE
    } else if(!missing(G)){
      Mat = G
      Zono =TRUE
    } else {
      print('No V-polytope or zonotope can be defined!')
      return(-1)
    }
    d = dim(Mat)[2] + 1
    m = dim(Mat)[1]
    b = rep(1, m)
    r = rep(0, d)
    r[1] = m
    r[2] = d
  } else {
    if (!missing(A)) {
      Mat = -A
      vec = b
      d = dim(Mat)[2] + 1
      m = dim(Mat)[1]
      r = rep(0,d)
      r[1] = m
      r[2] = d
    } else {
      print('matrix A is missing to define a H-polytope!')
      return(-1)
    }
  }
  Mat = matrix(cbind(b, Mat), ncol = dim(Mat)[2] + 1)
  Mat = matrix(rbind(r, Mat), ncol = dim(Mat)[2])
  
  # set flag for verbose mode
  verbose=FALSE
  
  # set flag for rotating only
  rotate_only = TRUE
  
  #---------------------#
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
  ball_only = FALSE
  #-------------------#
  
  Mat = vol_R(Mat, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only, 
              rotate_only, ball_only, sample_only, numpoints, variance, coordinate, rounding, verbose)
  
  # get elements "matrix" and "vector"
  # remove first row
  Mat = Mat[-c(1),]
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  if (vpoly || Zono){
    # in V-polytope or zonotope cases return only the marix
    return(A)
  } else {
    retList = list("A"=A, "b"=b)
    return(retList)
  }
}
