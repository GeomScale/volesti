#' Sample points from a convex Polytope (H-polytope, V-polytope or a zonotope) or use direct methods for uniform sampling from unit simplex and hypersphere
#'
#' Sample N points from a H or a V-polytope or a zonotope with uniform or spherical gaussian -centered in an internal point- target distribution.
#' #' The \eqn{d}-dimensional unit simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i\leq 1}, \eqn{x_i\geq 0}. The \eqn{d}-dimensional canonical simplex is the set of points \eqn{\vec{x}\in \R^d}, s.t.: \eqn{\sum_i x_i = 1}, \eqn{x_i\geq 0}.
#' 
#' @param P A convex polytope. It is an object from class (a) HPolytope or (b) VPolytope or (c) Zonotope.
#' @param N The number of points that the function is going to sample from the convex polytope. Default value is \eqn{100}.
#' @param distribution Optional. A list that contains parameters for the target distribution. Default distribution is uniform.
#' \itemize{
#'  \item{gaussian }{A boolean parameter. It declares spherical gaussian distribution as the target distribution. Default value is false.}
#'  \item{variance }{The variance for the spherical gaussian distribution. Default value is \eqn{1}.}
#' }
#' @param method Optional. A list that contains parameters for the random walk method. Default method is Coordinate Hit-and-Run.
#' \itemize{
#'  \item{direct }{A boolean parameter. It should be used for uniform sampling from the boundary or the interior of a hypersphere or from a unit or an arbitrary simplex. The arbitrary simplex has to be given as a V-polytope and the dimension should not be declared through method list. For the rest well known convex bodies it has to be declared the dimension and the type of body (simplex, sphere, ball).}
#'  \item{dim }{An integer that declares the dimension when direct flag is enabled for a unit simplex or a hypersphere (boundary or interior).}
#'  \item{body }{A string to request uniform sampling: (a) "simplex" to sample from an arbitrary simplex (when the simplex is given as a V-polytope) or a unit simplex (when no polytope is given and the dimension is declared), (b) "sphere" to sample from the boundary of a {d}-dimensional hypersphere centered at the origin and (c) to sample from the interior of the \eqn{d}-dimensional hypersphere centered at the origin. For (b) and (c) dimension should be given as well through method list.}
#'  \item{radius }{The radius of the \eqn{d}-dimensional hypersphere. Default value is \eqn{1}.}
#'  \item{WalkT }{A string to declare the random walk method: (a)"hnr" for Hit-and-Run or (b) "bw" for ball walk. Default method is Hit-and-Run.}
#'  \item{coord }{A boolean parameter for the hit-and-run. True for Coordinate Directions HnR, false for Random Directions HnR. Default value is TRUE.}
#'  \item{delta }{Optional. The radius for the ball walk.}
#'  \item{W }{Optional. The number of the steps for the random walk. Default value is \eqn{\lfloor 10+d/10\rfloor}.}
#' }
#' @param InnerPoint A \eqn{d}-dimensional numerical vector that defines a point in the interior of polytope P.
#' 
#' @references \cite{R.Y. Rubinstein and B. Melamed,
#' \dQuote{Modern simulation and modeling} \emph{ Wiley Series in Probability and Statistics,} 1998.}
#' @references \cite{A Smith, Noah and W Tromble, Roy,
#' \dQuote{Sampling Uniformly from the Unit Simplex,} \emph{ Center for Language and Speech Processing Johns Hopkins University,} 2004.}
#' @references \cite{Art B. Owen,
#' \dQuote{Monte Carlo theory, methods and examples,} \emph{ Copyright Art Owen,} 2009-2013.}
#' 
#' @return A \eqn{d\times N} matrix that contains, column-wise, the sampled points from the convex polytope.
#' @examples 
#' # uniform distribution from a 3d cube in V-representation using ball walk
#' P = GenCube(3, 'V')
#' points = sample_points(P, method = list("WalkT"="bw", "W"=5))
#' 
#' # gaussian distribution from a 2d unit simplex in H-representation with variance = 2
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' P = HPolytope(A=A, b=b)
#' points = sample_points(P, distribution = list("gaussian"=TRUE, "variance"=2))
#' @export
RRsample_points <- function(P, N, distribution, method, InnerPoint){
  
  if (!missing(P)) {
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
    walk_length = 10 + floor( dimension / 10 )
  } else {
    if(missing(method)) {
      stop("Neither a polytope or a direct method of sampling is declared.")
    }
    Mat = matrix(c(0,0))
    Zono = FALSE
    vpoly = FALSE
    walk_length = 0
  }
  
  gaussian = FALSE
  ball_walk = FALSE
  delta = -1
  coord = TRUE
  sam_simplex = FALSE
  sam_can_simplex = FALSE
  sam_arb_simplex = FALSE
  sam_ball = FALSE
  sam_sphere = FALSE
  numpoints = 100
  d = 0
  variance = 1.0
  if(!missing(method)) {
    if(!is.null(method$direct)) {
      if (method$direct) {
        if(!is.null(method$coord) || !is.null(method$W) || !is.null(method$WalkT) || !is.null(method$delta)) {
          warning("When direct sampling is enabled, the parameters of random walk should not be declared.")
        }
        if(is.null(method$body)) {
          stop("Direct sampling should be called only for simplices or hyperspheres.")
        }
        if (method$body=="simplex") {
          if(is.null(method$dim)){
            if(missing(P)) {
              stop("For direct sampling from a simplex you have to give the dimension of the unit simplex or a simplex in V-representation.")
            } else if(!vpoly) {
              stop("To sample uniformly points from a simplex, it has to be given in V-representation.")
            } else if(dimension != dim(Mat)[1]-2) {
              stop("The polytope is not a simplex.")
            }
            d = dimension
            sam_arb_simplex = TRUE
          } else {
            d = method$dim
            sam_simplex = TRUE
          }
        } else if(method$body=="can_simplex") {
          if(is.null(method$dim)) {
            stop("For direct sampling from a canonical simplex you have to give the dimension.")
          }
          d = method$dim
          sam_can_simplex = TRUE
        } else if (method$body=="sphere") {
          if(is.null(method$dim)) {
            stop("For direct sampling from a hypersphere you have to give the dimension of the hypersphere.")
          }
          d = method$dim
          sam_sphere = TRUE
        } else if(method$body=="ball") {
          if(is.null(method$dim)) {
            stop("For direct sampling from a ball you have to give the dimension of the ball.")
          }
          d = method$dim
          sam_ball = TRUE
        }
      }
    } else {
      if(missing(P)) {
        stop("No Polytope is defined.")
      }
      coordinate = TRUE
      if(!is.null(method$coord)){
        coordinate = method$coord
      }
      if(!is.null(method$WalkT)) {
        if(method$WalkT=="hnr") {
          ball_walk = FALSE
          delta = -1
        } else if(method$WalkT=="bw") {
          if(!is.null(method$coord)){
            warning("Ball walk and coordinate are both declared. Ball walk is going to be used.")
          }
          coordinate = TRUE
          ball_walk = TRUE
          delta = -1
          if(!is.null(method$delta)){
            delta = method$delta
          }
        } else {
          stop("Not a known random walk method.")
        }
      }
      if(!is.null(method$W)) {
        walk_length = method$W
        if (walk_length<=0) {
          stop("Walk length must be a positive value.")
        }
      }
    }
  }
  if(!missing(distribution)) {
    if (sam_simplex || sam_can_simplex || sam_arb_simplex || sam_ball || sam_sphere) {
      warning("The direct sampling can be used only for uniform points. Argument distribution should be NULL.")
    }
    if(!is.null(distribution$gaussian)) {
      gaussian = distribution$gaussian
      if(!is.null(distribution$uniform)) {
        if(gaussian && distribution$uniform) {
          stop("Only one target distribution has to be picked.")
        }
      }
      if(!is.null(distribution$variance)) {
        variance = distribution$variance
        if (variance<=0) {
          stop("Variance must be a positive value.")
        }
      }
    } else {
      if(!is.null(distribution$variance)) {
        warning("Variance should be declared when gaussian distribution is enabled.")
      }
    }
    
  }
  
  if(!missing(N)){
    numpoints = N
    if(numpoints<=0) {
      stop("Number of points has to be a positive integer.")
    }
  }
  
  # set too large vector for internal point if it is not given as an input
  innerpoint = rep(0,dim(Mat)[2] + 5)
  if (!missing(InnerPoint)) {
    if (sam_simplex || sam_can_simplex || sam_arb_simplex || sam_ball || sam_sphere) {
      warning("Direct sampling is enabled. Argument InnerPoint should be NULL.")
    }
    innerpoint = InnerPoint
  }
  
  points = Rsample_points(Mat, walk_length, innerpoint, gaussian,
                          ball_walk, delta, coord, vpoly, Zono,
                          sam_simplex, sam_can_simplex, sam_arb_simplex,
                          sam_ball, sam_sphere, numpoints, d, variance)
  
  return(points)
}
