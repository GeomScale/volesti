#' The main R function for volume approximation of a convex H or V Polytope
#'
#' For the volume approximation can be used two algorithms. Either volesti or CV. A H-polytope with m facets is described by a \eqn{m\times d} matrix A and a d-dimensional vector b, s.t.: \eqn{Ax\leq b}. A V-polytope is described as a set of d-dimensional points.
#'
#' @param list("argument"=value) A list that includes parameters for the chosen algorithm.
#' @param path The path to an ine or ext file that describes the H or V polytope respectively. If path is given then "matrix" and "vector" inputs are not needed.
#' @param matrix The matrix of the H polytope or the matrix that contains all the vertices of a V polytope row-wise. If the matrix is in ine file, for H-polytopes only (see examples), then the "vector" input is not needed.
#' @param vector Only for H-polytopes. The d-dimensional vector b that containes the constants of the facets.
#' @param walk_length Optional. The number of the steps for the random walk, default is \eqn{\lfloor 10 + d/10\rfloor}.
#' @param error Optional. Declare the goal for the approximation error. Default is 1 for volesti and \eqn{0.2} for CV.
#' @param Chebychev Optional. A d+1 vector that containes the chebychev center. The first d coordinates corresponds to the center and the last one to the radius of the chebychev ball.
#' @param CV Optional. A boolean parameter to use CV algorithm. Default value is false.
#' @param win_len Optional. The size of the window for the ratios' approximation in CV algorithm. Default value is \eqn{4 \ dimension^2 + 500}.
#' @param C Optional. a constant for the lower boud of \eqn{variance/mean^2} in schedule annealing.
#' @param N optional. The number of points we sample in each step of schedule annealing in CV algorithm. Default value is \eqn{500C + dimension^2 / 2}.
#' @param ratio Optional. parameter of schedule annealing, larger ratio means larger steps in schedule annealing. Default value is \eqn{1 - 1/dimension}.
#' @param frac Optional. the fraction of the total error to spend in the first gaussian. Default value is \eqn{0.1}.
#' @param ball_walk Optional. Boolean parameter to use ball walk, only for CV algorithm .Default value is false.
#' @param delta Optional. The radius for the ball walk.
#' @param verbose Optional. A boolean parameter for printing. Default is false.
#' @param vpoly A boolean parameter, has to be true when a V-polytope is given as input. Default value is false.
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is true.
#' @param rounding Optional. A boolean parameter to activate the rounding option. Default value is false.
#' 
#' @references \cite{I.Z.Emiris and V. Fisikopoulos,
#' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2014.}, 
#' @references \cite{B. Cousins and S. Vempala, \dQuote{A practical volume algorithm,} \emph{Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society,} 2015.}
#' 
#' @return The approximation of the volume of a convex H or V polytope.
#' @examples
#' # calling volesti algorithm for a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' vol = volume(list("matrix"=A, "vector"=b))
#' 
#' # calling CV algorithm for a V-polytope (3d cube)
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' vol = volume(list("matrix"=V, "CV"=TRUE, "Vpoly"=TRUE))
#' 
#' # a 2d unit simplex in H-representation using ine format matrix, calling volesti algorithm
#' A = matrix(c(3,3,0,0,-1,0,0,0,-1,1,1,1), ncol=3, nrow=4, byrow=TRUE)
#' vol = volume(list("matrix"=A))
volume <- function(Inputs){
  
  # flag for V-polytope
  Vpoly = FALSE
  if (!is.null(Inputs$Vpoly)) {
    Vpoly = Inputs$Vpoly
  }
  # initialization of the Polytope
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
      b = rep(1, m)
      r = rep(0, d)
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
  A = matrix(cbind(b,A), ncol = dim(A)[2] + 1)
  A = matrix(rbind(r,A), ncol = dim(A)[2])
  
  # set a too large vector for chebychev ball if it is not given as input
  Cheb_ball = rep(0, dim(A)[2] + 5)
  if (!is.null(Inputs$Chebychev)) {
    Cheb_ball = Inputs$Chebychev
  }
  
  # set flag for CV algorithm
  annealing = FALSE
  if (!is.null(Inputs$CV)) {
    annealing = Inputs$CV
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
  
  # set flag for rounding
  rounding = FALSE
  if (!is.null(Inputs$rounding)) {
    rounding = Inputs$rounding
  }
  
  # set the number of steps for the random walk
  if (!is.null(Inputs$walk_length)) {
    W=Inputs$walk_length
  } else {
    if (annealing) {
      W = 1
    }else{
      W = 10 + floor( dim(A)[2]/10 )
    }
  }
  
  # set the requested error
  if (!is.null(Inputs$error)) {
    e = Inputs$error
  } else {
    if (annealing) {
      e = 0.2
    } else {
      e = 1
    }
  }
  dimension = dim(A)[2] - 1
  
  # [CV] initialization
  win_len = 4 * ( dimension ^ 2 ) + 500
  if (!is.null(Inputs$win_len)) {
    win_len = Inputs$win_len
  }
  
  C = 2
  if (!is.null(Inputs$C)) {
    C = Inputs$C
  }
  ratio = 1 - 1 / dimension
  if (!is.null(Inputs$ratio)) {
    ratio = Inputs$ratio
  }
  N = 500 * C + ( dimension^2 ) / 2
  if (!is.null(Inputs$N)) {
    N = Inputs$N
  }
  frac = 0.1
  if (!is.null(Inputs$frac)) {
    frac = Inputs$frac
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

  #------------------------#
  round_only = FALSE
  rotate_only = FALSE
  sample_only = FALSE
  variance = 0
  numpoints = 0
  #-----------------------#
  # set the timer
  tim = proc.time()
  
  vol = vol_R(A, W, e, Cheb_ball, annealing, win_len, N, C, ratio, frac, ball_walk,
              delta, Vpoly, round_only, rotate_only, sample_only, numpoints, variance,
              coordinate, rounding, verbose)
  
  tim = proc.time()-tim
  if (verbose) {
    print(paste0('Total time: ', as.numeric(as.character(tim[3]))))
  }
  return(vol[1,1])
  
}
