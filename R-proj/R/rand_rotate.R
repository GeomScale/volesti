#' Apply a random rotation to a convex H or V-polytope.
#' 
#' Given a convex H or V polytope or a zonotope as input this function applies a random rotation.
#' 
#' @param list("argument"=value) A list that includes elements that describe the convex body that is given as input.
#' @param path The path to an ine (H-polytope) or ext (V-polytope, zonotope) file that describes the polytope. If path is given then "matrix" and "vector" inputs are not needed.
#' @param matrix The matrix of the H polytope or the matrix that containes all the d-dimensional vertices of a V polytope row-wise or a matrix that containes the d-dimensional segments that define a zonotope row-wise. If the matrix is in ine format, for H-polytopes only (see examples), then the "vector" input is not needed.
#' @param vector Only for H-polytopes. The m-dimensional vector b that containes the constants of the facets.
#' @param Vpoly A boolean parameter, has to be true when a V-polytope is given as input. Default value is false.
#' @param Zonotope A boolean parameter, has to be true when a zonotope is given as input. Default value is false.
#' @param verbose Optional. A boolean parameter for printing. Default is false.
#' 
#' @return A random rotation of the polytope that is given as an input. The output for a H-polytope is a list that containes elements "matrix" and "vector". For a V-polytope the output is a \eqn{m\times d} matrix that contains the \eqn{m} d-dimensional vertices of the V polytope row-wise. For a zonotope is a \eqn{m\times d} matrix that containes the \eqn{m} d-dimensional segments row-wise.
#' @examples
#' # rotate a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' listHpoly = rand_rotate(list("matrix"=A, "vector"=b))
#' 
#' # rotate a V-polytope (3d cube)
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' matVpoly = rand_rotate(list("matrix"=V, "Vpoly"=TRUE))
#' 
#' # rotate a 5-dimensional zonotope defined by the Minkowski sum of 15 segments
#' Zono = GenZonotope(5,15)
#' MatZono = rand_rotate(list("matrix"=Zono, "Zonotope"=TRUE))
rand_rotate <- function(Inputs){
  
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
  
  # set flag for verbose mode
  verbose=FALSE
  if(!is.null(Inputs$verbose)){
    verbose=Inputs$verbose
  }
  
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
  #-------------------#
  
  Mat = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              Vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen, round_only, 
              rotate_only, sample_only, numpoints, variance, coordinate, rounding, verbose)
  
  # get elements "matrix" and "vector"
  retList = modifyMat(Mat)
  if (Vpoly){
    # in V-polytope case return only the marix
    return(retList$matrix)
  } else {
    return(retList)
  }
}
