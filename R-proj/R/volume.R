#' The main R function for volume approximation of a convex Polytope (H-polytope, V-polytope or a zonotope)
#'
#' For the volume approximation can be used two algorithms. Either SequenceOfBalls or CoolingGaussian. A H-polytope with \eqn{m} facets is described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{Ax\leq b}. A V-polytope is described as a set of \eqn{d}-dimensional points. A zonotope is desrcibed by the Minkowski sum of \eqn{d}-dimensional segments.
#'
#' @param A Only for H-polytopes. The \eqn{m\times d} matrix \eqn{A} that containes the directions of the \eqn{m} facets.
#' @param b Only for H-polytopes. The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets s.t.: \eqn{Ax\leq b}.
#' @param V Only for V-polytopes. The \eqn{m\times d} matrix V that containes row-wise the \eqn{m} \eqn{d}-dimensional vertices of the polytope.
#' @param G Only for zonotopes. The \eqn{m\times d} matrix G that containes row-wise the \eqn{m} \eqn{d}-dimensional segments that define a zonotope.
#' @param walk_length Optional. The number of the steps for the random walk. Default value is \eqn{\lfloor 10 + d/10\rfloor}.
#' @param error Optional. Declare the goal for the approximation error. Default value is \eqn{1} for SequenceOfBalls and \eqn{0.2} for CoolingGaussian.
#' @param InnerVec Optional. A \eqn{d+1} vector that containes an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an internal ball.
#' @param CG Optional. A boolean parameter to use CoolingGaussian algorithm. Default value is false.
#' @param win_len Optional. The size of the window for the ratios' approximation in CG algorithm. Default value is \eqn{4 \cdot dimension^2 + 500}.
#' @param C Optional. A constant for the lower bound of \eqn{variance/mean^2} in schedule annealing of CG algorithm.
#' @param N optional. The number of points we sample in each step of schedule annealing in CG algorithm. Default value is \eqn{500C + dimension^2 / 2}.
#' @param ratio Optional. Parameter of schedule annealing of CG algorithm, larger ratio means larger steps in schedule annealing. Default value is \eqn{1 - 1/dimension}.
#' @param frac Optional. The fraction of the total error to spend in the first gaussian in CG algorithm. Default value is \eqn{0.1}.
#' @param ball_walk Optional. Boolean parameter to use ball walk. Default value is false.
#' @param delta Optional. The radius for the ball walk.
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is true.
#' @param rounding Optional. A boolean parameter to activate the rounding option. Default value is false.
#' 
#' @references \cite{I.Z.Emiris and V. Fisikopoulos,
#' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2014.}, 
#' @references \cite{B. Cousins and S. Vempala, \dQuote{A practical volume algorithm,} \emph{Springer-Verlag Berlin Heidelberg and The Mathematical Programming Society,} 2015.}
#' 
#' 
#' @return The approximation of the volume of a convex polytope.
#' @examples
#' # calling SOB algorithm for a H-polytope (2d unit simplex)
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' vol = volume(A=A, b=b)
#' 
#' # calling CG algorithm for a V-polytope (3d cube)
#' Vmat = GenSimplex(2,'V')
#' vol = volume(V=Vmat, CG=TRUE)
#' 
#' # calling CG algorithm for a 5-dimensional zonotope defined as the Minkowski sum of 10 segments
#' zonotope = GenZonotope(2, 4)
#' vol = volume(G=zonotope, rounding=TRUE, CG=TRUE)
#' @export
#' @useDynLib volesti
#' @importFrom Rcpp evalCpp
#' @importFrom "utils" "read.csv"
#' @exportPattern "^[[:alpha:]]+"
volume <- function(P, walk_length, error, InnerBall, Algo, WalkType, rounding){
  
  repr = class(P)[1]
  if (repr == "HPolytope") {
    vpoly = FALSE
    Zono = FALSE
  } else if(repr == "VPolytope") {
    vpoly = TRUE
    Zono = FALSE
  } else if(repr == "Zonotope") {
    vpoly = FALSE
    Zono = TRUE
  } else {
    stop("Not a known polytope representation.")
  }
  
  Mat = P$get_mat()
  
  dimension = dim(Mat)[2] - 1
  
  CG = FALSE
  win_len = 0
  C = 0
  N = 0
  ratio = 0
  frac = 0
  if (!missing(Algo)) {
    if(!is.null(Algo$CG)) {
      if(Algo$CG) {
        CG = TRUE
        
        win_len=4*(dimension^2)+500
        if(!is.null(Algo$window_len)){
          win_len=Algo$window_len
        }
        C=2
        if(!is.null(Algo$C)){
          C=Algo$C
        }
        ratio=1-1/dimension
        if(!is.null(Algo$ratio)){
          ratio=Algo$ratio
        }
        N=500*C+(dimension^2)/2
        if(!is.null(Algo$N)){
          N=Algo$N
        }
        frac=0.1
        if(!is.null(Algo$frac)){
          frac=Algo$frac
        }
      } else if (!is.null(Algo$SOB)) {
        if (!Algo$SOB) {
          stop("You have to choose between two Algorithms!")
        }
      } else {
        warning("CG is false and no flag for SOB. The latest algorithm will be used.")
      }
    }
  }
  
  # set the number of steps for the random walk
  if (!missing(walk_length)) {
    W = walk_length
  } else {
    if (CG) {
      W = 1
    }else{
      W = 10 + floor( dimension / 10 )
    }
  }
  
  # set the requested error
  if (!missing(error)) {
    e = error
  } else {
    if (CG) {
      e = 0.2
    } else {
      e = 1
    }
  }
  
  # set a too large vector for chebychev ball if it is not given as input
  InnerB = rep(0, dimension + 5)
  if (!missing(InnerBall)) {
    InnerB = InnerBall
  }
  
  round = FALSE
  if (!missing(rounding)) {
    round = rounding
  }
  
  if (missing(WalkType)) {
    ball_walk = FALSE
    delta = -1
    coordinate = TRUE
  } else {
    if(is.null(WalkType$method)){
      stop("No method for random wak was picked.")
    }
    if(WalkType$method=="hnr") {
      ball_walk = FALSE
      delta = -1
      coordinate = TRUE
      if(!is.null(WalkType$coordinate)){
        coordinate = WalkType$coordinate
      }
    } else if(WalkType$method=="bw") {
      coordinate = TRUE
      ball_walk = TRUE
      delta = -1
      if(!is.null(WalkType$delta)){
        delta = WalkType$delta
      }
    } else {
      stop("Not a known random walk method.")
    }
  }
  
  vol = Rvolume(Mat, W, e, InnerB, CG, win_len, N, C, ratio,
                frac, ball_walk, delta, vpoly, Zono, coordinate, round)

  return(vol)
  
}
