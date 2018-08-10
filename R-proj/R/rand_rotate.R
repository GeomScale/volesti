#' Apply a random rotation to a convex H or V-polytope.
#' 
#' Given a convex H or V polytope as input and then a random rotation is applied to the polytope.
#' 
#' @param list("argument"=value) A list that includes elements that describe the convex polytope given as input.
#' @param path The path to an ine or ext file that describes the H or V polytope respectively. If path is given then "matrix" and "vector" inputs are not needed.
#' @param matrix The matrix of the H polytope or the matrix that contains all the vertices of a V polytope row-wise. If the matrix is in ine file, for H-polytopes only (see examples), then the "vector" input is not needed.
#' @param vector Only for H-polytopes. The d-dimensional vector b that containes the constants of the facets.
#' @param vpoly A boolean parameter, has to be true when a V-polytope is given as input. Default value is false.
#' @param verbose Optional. A boolean parameter for printing. Default is false.
#' 
#' @return A H or V-polytope which is a random rotation of the polytope that is given as an input. The output for a H-polytope is a list that containes elements "matrix" and "vector". For a V-polytope the output is a \eqn{k\times d} matrix that contains the \eqn{k} vertices of the V polytope row-wise.
#' @examples
#' #rotate a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' listHpoly = rand_rotate(list("matrix"=A, "vector"=b))
#' 
#' #rotate a V-polytope (3d cube)
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' matVpoly = rand_rotate(list("matrix"=V, "Vpoly"=TRUE))
rand_rotate <- function(Inputs){
  
  # set flag for V-polytope
  Vpoly = FALSE
  if (!is.null(Inputs$Vpoly)) {
    Vpoly = Inputs$Vpoly
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
    if (Vpoly) {
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
    if (Vpoly) {
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
  #-------------------#
  
  Mat = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk, delta,
              Vpoly, round_only, rotate_only, sample_only, numpoints, variance, coordinate,
              rounding, verbose)
  
  # get elements "matrix" and "vector"
  retList = modifyMat(Mat)
  if (Vpoly){
    # in V-polytope case return only the marix
    return(retList$matrix)
  } else {
    return(retList)
  }
}
