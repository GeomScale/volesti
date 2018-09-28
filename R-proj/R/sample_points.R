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
sample_points <- function(P, N, distribution, method, InnerPoint){
  
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
            } else if(!vpoly || Zono) {
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
