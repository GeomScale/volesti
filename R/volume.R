#' The main R function for volume approximation of a convex Polytope (H-polytope, V-polytope or a zonotope).
#'
#' For the volume approximation can be used two algorithms. Either SequenceOfBalls or CoolingGaussian. A H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}. A V-polytope is described as a set of \eqn{d}-dimensional points. A zonotope is desrcibed by the Minkowski sum of \eqn{d}-dimensional segments.
#'
#' @param A Only for H-polytopes. The \eqn{m\times d} matrix \eqn{A} that containes the directions of the \eqn{m} facets.
#' @param b Only for H-polytopes. The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets s.t.: \eqn{Ax\leq b}.
#' @param V Only for V-polytopes. The \eqn{m\times d} matrix V that containes row-wise the \eqn{m} \eqn{d}-dimensional vertices of the polytope.
#' @param G Only for zonotopes. The \eqn{m\times d} matrix G that containes row-wise the \eqn{m} \eqn{d}-dimensional segments that define a zonotope.
#' @param walk_length Optional. The number of the steps for the random walk, default is \eqn{\lfloor 10 + d/10\rfloor}.
#' @param error Optional. Declare the goal for the approximation error. Default is \eqn{1} for SequenceOfBalls and \eqn{0.2} for CoolingGaussian.
#' @param InnerVec Optional. A \eqn{d+1} vector that containes an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined as the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,m}.
#' @param CG Optional. A boolean parameter to use CoolingGaussian algorithm. Default value is false.
#' @param win_len Optional. The size of the window for the ratios' approximation in CG algorithm. Default value is \eqn{4 \ dimension^2 + 500}.
#' @param C Optional. A constant for the lower bound of \eqn{variance/mean^2} in schedule annealing of CG algorithm.
#' @param N optional. The number of points we sample in each step of schedule annealing in CG algorithm. Default value is \eqn{500C + dimension^2 / 2}.
#' @param ratio Optional. Parameter of schedule annealing of CG algorithm, larger ratio means larger steps in schedule annealing. Default value is \eqn{1 - 1/dimension}.
#' @param frac Optional. The fraction of the total error to spend in the first gaussian in CG algorithm. Default value is \eqn{0.1}.
#' @param ball_walk Optional. Boolean parameter to use ball walk. Default value is false.
#' @param delta Optional. The radius for the ball walk.
#' @param verbose Optional. A boolean parameter for printing. Default is false.
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
#' vol = volume(A=A, b=b)
#' 
#' # calling CV algorithm for a V-polytope (3d cube)
#' V = matrix(c(-1,1,-1,-1,-1,1,-1,1,1,-1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,-1,-1), ncol=3, nrow=8, byrow=TRUE)
#' vol = volume(V=V, CG=TRUE)
#' 
#' # calling Gaussian-Cooling algorithm for a 5-dimensional zonotope defined as the Minkowski sum of 10 segments
#' zonotope = GenZonotope(5, 10)
#' vol = volume(G=zonotope, rounding=TRUE, CG=TRUE)
volume <- function(A, b, V, G, walk_length, error, InnerVec, CG, win_len,
                   C, N, ratio, frac, ball_walk, delta, verbose, 
                   coordinate, rounding) {
  
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
  
  dimension = dim(Mat)[2] - 1
  
  # set a too large vector for chebychev ball if it is not given as input
  InnerBall = rep(0, dimension + 5)
  if (!missing(InnerVec)) {
    InnerBall = InnerVec
  }
  
  # set flag for CV algorithm
  annealing = FALSE
  if (!missing(CG)) {
    annealing = CG
  }
  
  # set flag for verbose mode
  verb = FALSE
  if (!missing(verbose)) {
    verb = verbose
  }
  
  # set flag for Coordinate or Random Directions HnR
  coord = TRUE
  if (!missing(coordinate)) {
    coord = coordinate
  }
  
  # set flag for rounding
  round = FALSE
  if (!missing(rounding)) {
    round = rounding
  }
  
  # set the number of steps for the random walk
  if (!missing(walk_length)) {
    W = walk_length
  } else {
    if (annealing) {
      W = 1
    }else{
      W = 10 + floor( dimension / 10 )
    }
  }
  
  # set the requested error
  if (!missing(error)) {
    e = error
  } else {
    if (annealing) {
      e = 0.2
    } else {
      e = 1
    }
  }
  
  
  # [CG] initialization
  window_len = 4 * ( dimension ^ 2 ) + 500
  if (!missing(win_len)) {
    window_len = win_len
  }
  
  c = 2
  if (!missing(C)) {
    c = C
  }
  Ratio = 1 - 1 / dimension
  if (!missing(ratio)) {
    Ratio = ratio
  }
  NN = 500 * c + ( dimension^2 ) / 2
  if (!missing(N)) {
    NN = N
  }
  Frac = 0.1
  if (!missing(frac)) {
    Frac = frac
  }
  
  # set flag for the ball walk
  ballwalk = FALSE
  if (!missing(ball_walk)) {
    ballwalk = ball_walk
  }
  
  # set the radius for the ball walk. Negative value means that is not given as input
  Delta = -1
  if (!missing(delta)) {
    Delta = delta
  }

  #------------------------#
  round_only = FALSE
  rotate_only = FALSE
  sample_only = FALSE
  variance = 0
  numpoints = 0
  gen_only = FALSE
  Vpoly_gen = FALSE
  kind_gen = -1
  dim_gen = 0
  m_gen = 0
  exact_zono = FALSE
  #-----------------------#
  # set the timer
  tim = proc.time()
  
  vol = vol_R(Mat, W, e, InnerBall, annealing, window_len, NN, c, Ratio, Frac, ballwalk,
              Delta, vpoly, Zono, exact_zono, gen_only, Vpoly_gen, kind_gen, dim_gen, m_gen,
              round_only, rotate_only, sample_only, numpoints, variance, coord,
              round, verb)
  
  tim = proc.time()-tim
  if (verb) {
    print(paste0('Total time: ', as.numeric(as.character(tim[3]))))
  }
  return(vol[1,1])
  
}
