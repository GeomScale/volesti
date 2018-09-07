#' Apply rounding to a convex polytope (H-polytope, V-polytope or a zonotope).
#' 
#' Given a convex H or V polytope or a zonotope as input this function computes a rounding based on minimum volume enclosing ellipsoid of a pointset.
#' 
#' @param A Only for H-polytopes. The \eqn{m\times d} matrix \eqn{A} that containes the directions of the \eqn{m} facets.
#' @param b Only for H-polytopes. The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets s.t.: \eqn{Ax\leq b}.
#' @param V Only for V-polytopes. The \eqn{m\times d} matrix V that containes row-wise the \eqn{m} \eqn{d}-dimensional vertices of the polytope.
#' @param G Only for zonotopes. The \eqn{m\times d} matrix G that containes row-wise the \eqn{m} \eqn{d}-dimensional segments that define a zonotope.
#' @param walk_length Optional. The number of the steps for the random walk, default is \eqn{\lfloor 10+d/10\rfloor}.
#' @param ball_walk Optional. Boolean parameter to use ball walk, only for CG algorithm. Default value is false.
#' @param delta Optional. The radius for the ball walk.
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is true.
#' @param verbose Optional. A boolean parameter for printing. Default is false.
#' 
#' @return Is a list that containes elements to describe the rounded polytope, i.e. "matrix" and "vector" for H-polytopes and just "matrix" for V-polytopes and zonotopes, containing the verices or segments row-wise. For both representations the list containes element "round_value" which is the determinant of the square matrix of the linear transformation that was applied on the polytope that is given as input.
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' listHpoly = round_polytope(A=A, b=b)
#' 
#' # rotate a V-polytope (3d cube) using Random Directions HnR
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' ListVpoly = round_polytope(V=V, coordinate=FALSE)
#' 
#' # rotate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' Zono = GenZonotope(10,20)
#' ListZono = round_polytope(G=Zono)
round_polytope <- function(A, b, V, G, walk_length, ball_walk, delta, coordinate, verbose) {
  
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
  
  # set the number of steps for the random walk
  W = 10 + floor((dim(Mat)[2] - 1) / 10)
  if (!missing(walk_length)) {
    W = walk_length
  }
  
  # set flag for Coordinate or Random Directions HnR
  coord = TRUE
  if (!missing(coordinate)) {
    coord = coordinate
  }
  
  # set flag for ball walk
  ballwalk = FALSE
  if (!missing(ball_walk)) {
    ballwalk = ball_walk
  }
  
  # set the radius for the ball walk. Negative value means that is not given as input
  Delta = -1
  if (!missing(delta)) {
    Delta = delta
  }
  
  # set flag for verbose mode
  verb = FALSE
  if (!missing(verbose)) {
    verb = verbose
  }
  
  #set round_only flag
  round_only = TRUE
  
  #---------------------#
  rotate_only = FALSE
  e = 0
  InnerBall = rep(0,dim(Mat)[2] + 5)
  annealing = FALSE
  win_len = 0
  N = 0
  C = 0
  ratio = 0
  frac = 0
  sample_only = FALSE
  numpoints = 0
  variance = 0
  rounding = FALSE
  gen_only = FALSE
  Vpoly_gen = FALSE
  kind_gen = -1
  dim_gen = 0
  m_gen = 0
  exact_zono = FALSE
  #---------------------#
  
  Mat = vol_R(Mat, W, e, InnerBall, annealing, win_len, N, C, ratio, frac, ballwalk,
              Delta, vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only,
              rotate_only, sample_only, numpoints, variance, coord, rounding, verb)
  # get first row which has the info for round_value
  r = Mat[c(1),]
  round_value = r[1]
  # remove first row
  Mat = Mat[-c(1),]
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  if (vpoly) {
    retList = list("V"=A, "round_value"=round_value)
    return(retList)
  }else if (Zono) {
    retList = list("G"=A, "round_value"=round_value)
    return(retList)
  } else {
    retList = list("A"=A, "b"=b, "round_value"=round_value)
    return(retList)
  }
}
