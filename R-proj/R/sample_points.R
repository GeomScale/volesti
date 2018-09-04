#' Sample points from a convex Polytope
#'
#' Sample N points from a H or a V-polytope with uniform or spherical gaussian target distribution.
#' @param list("argument"=value) A list that includes parameters for the chosen target distribution and the random walk algorithm.
#' @param path The path to an ine or ext file that describes the H or V polytope respectively. If path is given then "matrix" and "vector" inputs are not needed.
#' @param matrix The matrix A of a H-polytope or the matrix V that contains all the vertices of a V polytope row-wise. If it is in ine format, only for H-polytopes, then the input "vector" is not needed.
#' @param vector Only for H-polytopes. The d-dimensional vector b that containes the constants of the facets.
#' @param walk_length Optional. The number of the steps for the random walk, default is \eqn{\lfloor 10+d/10\rfloor}.
#' @param internal_point Optional. A d-dimensional vector that containes an internal point of the polytope.
#' @param gaussian Optional. A boolean parameter to sample with gaussian target distribution. Default value is false.
#' @param variance Optional. The variance for the spherical gaussian. Default value is \eqn{1}.
#' @param N The number of points that the function is going to sample from the convex polytope. Default value is \eqn{100}.
#' @param ball_walk Optional. Boolean parameter to use ball walk for the sampling. Default value is false.
#' @param delta Optional. The radius for the ball walk.
#' @param verbose Optional. A boolean parameter for printing. Default is false.
#' @param Vpoly A boolean parameter, has to be true when a V-polytope is given as input. Default value is false.
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is true.
#' @return A \eqn{d\times N} matrix that contains, column-wise, the sampled points from the convex polytope.
#' @examples 
#' #uniform distribution from a 3d cube described by a set of vertices
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' points = sample_points(list("matrix"=V, "Vpoly"=TRUE, "N"=1000))
#' 
#' #gaussian distribution from a 2d unit simplex in H-representation with variance = 2
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' points = sample_points(list("matrix"=A, "vector"=b, "gaussian"=TRUE, "variance"=2))
sample_points <- function(Inputs){
  
  #set flag for V-polytope
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
  
  # set the number of points to sample
  N = 100
  if (!is.null(Inputs$N)) {
    N = Inputs$N
  }
  
  # set too large vector for internal point if it is not given as an input
  internal_point = rep(0,dim(A)[2] + 5)
  if (!is.null(Inputs$internal_point)) {
    internal_point = Inputs$internal_point
  }
  
  # set flag for spherical gaussian distribution
  gaussian = FALSE
  if (!is.null(Inputs$gaussian)) {
    gaussian = Inputs$gaussian
  }
  
  # set variance
  variance = 1
  if (!is.null(Inputs$variance)) {
    variance = Inputs$variance
  }
  
  # set flag for verbose mode
  verbose = FALSE
  if (!is.null(Inputs$verbose)) {
    verbose = Inputs$verbose
  }
  
  # set flag for Coordinate or Random Directions HnR
  coordinate = TRUE
  if (!is.null(Inputs$coordinate)) {
    coordinate = Inputs$coordinate
  }
  
  # set the number of steps for the random walk
  W = 10 + floor((dim(A)[2] - 1) / 10)
  if (!is.null(Inputs$walk_length)) {
    W = Inputs$walk_length
  }
  
  # set flag for the ball walk
  ball_walk = FALSE
  if (!is.null(Inputs$ball_walk)) {
    ball_walk = Inputs$ball_walk
  }
  
  # set the radius for the ball walk. Negative value means that is not given as input
  delta = -1
  if (!is.null(Inputs$delta)) {
    delta = Inputs$delta
  }
  
  #--------------------#
  rounding = FALSE
  e = 0
  win_len = 0
  C = 0
  ratio = 0
  NN = 0
  frac = 0
  
  round_only = FALSE
  rotate_only = FALSE
  sample_only = TRUE
  gen_only = FALSE
  Vpoly_gen = FALSE
  kind_gen = -1
  dim_gen = 0
  m_gen = 0
  exact_zono = FALSE
  #---------------------#
  
  # set timer
  tim = proc.time()
  
  points = vol_R(A, W, e, internal_point, gaussian, win_len, NN, C, ratio, frac,
                 ball_walk, delta, Vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen,
                 m_gen, round_only, rotate_only, sample_only, N, variance, coordinate,
                 rounding, verbose)
  
  tim = proc.time() - tim
  if (verbose) {
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  return(points)
}
