#' Apply rounding to a convex polytope (H-polytope, V-polytope or a zonotope).
#' 
#' Given a convex H or V polytope or a zonotope as input this function computes a rounding based on minimum volume enclosing ellipsoid of a pointset.
#' 
#' @param list("argument"=value) A list that includes parameters for the rounding.
#' @param path The path to an ine (H-polytope) or ext (V-polytope, zonotope) file that describes the polytope. If path is given then "matrix" and "vector" inputs are not needed.
#' @param matrix The \eqn{m\times d} matrix \eqn{A} of the H polytope or the \eqn{m\times d} matrix that containes all the \eqn{m} \eqn{d}-dimensional vertices of a V-polytope row-wise or a \eqn{m\times d} matrix that containes all the \eqn{m} \eqn{d}-dimensional segments that define a zonotope row-wise. If the matrix is in ine format, for H-polytopes only (see \eqn{volume} function example), then the "vector" input is not needed.
#' @param vector Only for H-polytopes. The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets, s.t.: \eqn{Ax\leq b}.
#' @param Vpoly A boolean parameter, has to be true when a V-polytope is given as input. Default value is false.
#' @param Zonotope A boolean parameter, has to be true when a zonotope is given as input. Default value is false.
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
#' listHpoly = round_polytope(list("matrix"=A, "vector"=b))
#' 
#' # rotate a V-polytope (3d cube) using Random Directions HnR
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' ListVpoly = round_polytope(list("matrix"=V, "Vpoly"=TRUE, "coordinate"=FALSE))
#' 
#' # rotate a 10-dimensional zonotope defined by the Minkowski sum of 20 segments
#' Zono = GenZonotope(10,20)
#' ListZono = round_polytope(list("matrix"=Zono, "Zonotope"=TRUE))
round_polytope <- function(Inputs){
  
  # set flag for V-polytope
  Vpoly = FALSE
  if (!is.null(Inputs$Vpoly)) {
    Vpoly = Inputs$Vpoly
  }
  Zono = FALSE
  if (!is.null(Inputs$Zonotope)) {
    Zono = Inputs$Zonotope
  }
  
  # polytope initialization
  if (!is.null(Inputs$path)) {
    A = ineToMatrix(read.csv(Inputs$path))
    r = A[1,]
    x = modifyMat(A)
    A = x$matrix
    b = x$vector
  } else if (!is.null(Inputs$vector)) {
    b = Inputs$vector
    A = -Inputs$matrix
    d = dim(A)[2] + 1
    m = dim(A)[1]
    r = rep(0,d)
    r[1] = m
    r[2] = d
  } else if (!is.null(Inputs$matrix)) {
    if (Vpoly || Zono) {
      A = Inputs$matrix
      d = dim(A)[2] + 1
      m = dim(A)[1]
      b = rep(1,m)
      r = rep(0,d)
      r[1] = m
      r[2] = d
    } else {
      r = Inputs$matrix[1,]
      x = modifyMat(Inputs$matrix)
      A = x$matrix
      b = x$vector
    }
  } else {
    if (Zono) {
      print('No Zonotope defined from input!')
    } else if (Vpoly) {
      print('No V-polytope defined from input!')
    } else {
      print('No H-polytope defined from input!')
    }
    return(-1)
  }
  A = matrix(cbind(b,A), ncol=dim(A)[2] + 1)
  A = matrix(rbind(r,A), ncol=dim(A)[2])
  
  # set the number of steps for the random walk
  W = 10 + floor((dim(A)[2] - 1) / 10)
  if (!is.null(Inputs$walk_length)) {
    W = Inputs$walk_length
  }
  
  # set flag for Coordinate or Random Directions HnR
  coordinate = TRUE
  if (!is.null(Inputs$coordinate)) {
    coordinate = Inputs$coordinate
  }
  
  # set flag for ball walk
  ball_walk = FALSE
  if (!is.null(Inputs$ball_walk)) {
    ball_walk = Inputs$ball_walk
  }
  
  # set the radius for the ball walk. Negative value means that is not given as input
  delta = -1
  if (!is.null(Inputs$delta)) {
    delta = Inputs$delta
  }
  
  # set flag for verbose mode
  verbose = FALSE
  if (!is.null(Inputs$verbose)) {
    verbose = Inputs$verbose
  }
  
  #set round_only flag
  round_only = TRUE
  
  #---------------------#
  rotate_only = FALSE
  e = 0
  Cheb_ball = rep(0,dim(A)[2] + 5)
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
  
  Mat = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk,
              delta, Vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only,
              rotate_only, sample_only, numpoints, variance, coordinate, rounding, verbose)
  # get first row which has the info for round_value
  r = Mat[c(1),]
  round_value = r[1]
  # get "matrix" and "vector" elements
  retList = modifyMat(Mat)
  if (Vpoly) {
    output = list("matrix"=retList$matrix, "round_value"=round_value)
    return(output)
  } else {
    output = list("matrix"=retList$matrix, "vector"=retList$vector, "round_value"=round_value)
    return(output)
  }
}
