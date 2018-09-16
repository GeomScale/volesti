#' Sample points from a convex Polytope (H-polytope, V-polytope or a zonotope)
#'
#' Sample N points from a H or a V-polytope or a zonotope with uniform or spherical gaussian -centered in an internal point- target distribution.
#' 
#' @param A Only for H-polytopes. The \eqn{m\times d} matrix \eqn{A} that containes the directions of the \eqn{m} facets.
#' @param b Only for H-polytopes. The \eqn{m}-dimensional vector \eqn{b} that containes the constants of the \eqn{m} facets s.t.: \eqn{Ax\leq b}.
#' @param V Only for V-polytopes. The \eqn{m\times d} matrix V that containes row-wise the \eqn{m} \eqn{d}-dimensional vertices of the polytope.
#' @param G Only for zonotopes. The \eqn{m\times d} matrix G that containes row-wise the \eqn{m} \eqn{d}-dimensional segments that define a zonotope.
#' @param walk_length Optional. The number of the steps for the random walk. Default value is \eqn{\lfloor 10+d/10\rfloor}.
#' @param internal_point Optional. A \eqn{d}-dimensional vector that containes the coordinates of an internal point of the polytope. If it is not given then for H-polytopes the Chebychev center is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev center of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we use the origin.
#' @param gaussian Optional. A boolean parameter to sample with gaussian target distribution. Default value is false.
#' @param variance Optional. The variance for the spherical gaussian. Default value is \eqn{1}.
#' @param N The number of points that the function is going to sample from the convex polytope. Default value is \eqn{100}.
#' @param ball_walk Optional. Boolean parameter to use ball walk for the sampling. Default value is false.
#' @param delta Optional. The radius for the ball walk.
#' @param verbose Optional. A boolean parameter for printing. Default value is false.
#' @param coordinate Optional. A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is true.
#' @return A \eqn{d\times N} matrix that containes, column-wise, the sampled points from the convex polytope.
#' @examples 
#' # uniform distribution from a 3d cube described by a set of vertices
#' Vmat = GenCube(3, 'V')
#' points = sample_points(V=Vmat, N=1000)
#' 
#' # gaussian distribution from a 2d unit simplex in H-representation with variance = 2
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' points = sample_points(A=A, b=b, gaussian=TRUE, variance=2)
#' @export
sample_points <- function(A, b, V, G, walk_length, internal_point,
                          gaussian, variance, N, ball_walk, delta,
                          verbose, coordinate){
  
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
  
  # set the number of points to sample
  n = 100
  if (!missing(N)) {
    n = N
  }
  
  # set too large vector for internal point if it is not given as an input
  internalpoint = rep(0,dim(Mat)[2] + 5)
  if (!missing(internal_point)) {
    internalpoint = internal_point
  }
  
  # set flag for spherical gaussian distribution
  Gaussian = FALSE
  if (!missing(gaussian)) {
    Gaussian = gaussian
  }
  
  # set variance
  var = 1
  if (!missing(variance)) {
    var = variance
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
  
  # set the number of steps for the random walk
  W = 10 + floor((dim(Mat)[2] - 1) / 10)
  if (!missing(walk_length)) {
    W = walk_length
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
  ball_only = FALSE
  sam_simplex = FALSE
  sam_can_simplex = FALSE
  sam_arb_simplex = FALSE
  #---------------------#
  
  # set timer
  tim = proc.time()
  
  points = vol_R(Mat, W, e, internalpoint, Gaussian, win_len, NN, C, ratio, frac,
                 ballwalk, Delta, vpoly, Zono, exact_zono, gen_only, Vpoly_gen,
                 kind_gen, dim_gen, m_gen, round_only, rotate_only, ball_only,
                 sample_only, sam_simplex, sam_can_simplex, sam_arb_simplex, n,
                 var, coord, rounding, verb)
  
  tim = proc.time() - tim
  if (verb) {
    print(paste0('Total time: ',as.numeric(as.character(tim[3]))))
  }
  return(points)
}
