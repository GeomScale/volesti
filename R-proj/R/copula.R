#' Construct a copula using uniform sampling from the unit simplex
#' 
#' Given two families of parallel hyperplanes (or a family of parallel hyperplanes and a family of concentric ellispoids centered at the origin) intersecting the canonical simplex, this function samples from the canonical simplex and construct an approximation of the bivariate probability distribution, called copula.
#' 
#' @param h1 A \eqn{d}-dimensional vector that describes the direction of the first family of parallel hyperplanes.
#' @param h2 A \eqn{d}-dimensional vector that describes the direction of the second family of parallel hyperplanes.
#' @param E The \eqn{d\times d} symmetric positive define matrix of the family of concentric ellipsoids centered at the origin.
#' @param numSlices The number of the slices for the copula. Default value is 100.
#' @param N The number of points to sample. Default value is \eqn{4\cdot 10^6}.
#'
#' @references \cite{L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos,
#' \dQuote{Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises,} \emph{Proc. of Symposium on Computational Geometry, Budapest, Hungary,} 2018.}
#'
#' @return A \eqn{numSlices\times numSlices} copula.
#' @examples 
#' # compute a copula for two families of parallel hyperplanes
#' h1 = runif(n = 10, min = 1, max = 1000)
#' h1 = h1 / 1000
#' h2=runif(n = 10, min = 1, max = 1000)
#' h2 = h2 / 1000
#' cop = copula(h1=h1, h2=h2, numSlices = 10, N = 100000)
#' 
#' # compute a copula for a family of parallel hyperplanes and a family of conentric ellipsoids
#' h1 = runif(n = 10, min = 1, max = 1000)
#' h1 = h1 / 1000
#' E = replicate(10, rnorm(20))
#' E = cov(E)
#' cop = copula(h1=h1, E=E, numSlices=10, N=100000)

copula <-function(h1, h2, E, numSlices, N) {
  if(missing(h1)) {
    print('Wrong inputs..see the documentaion')
    return(0)
  } else if (missing(E)) {
    if (missing(h2)) {
      print('Wrong inputs..see the documentaion')
      return(0)
    }
    if (length(h1)!=length(h2)  || length(h1)<2){
      print('Wrong inputs..see the documentaion')
      return(0)
    }
    E=matrix(c(0,0))
  } else {
    if (!missing(h2)) {
      print('Wrong inputs..see the documentaion')
      return(0)
    }
    if(length(h1)!=dim(E)[2] || length(h1)<2) {
      print('Wrong inputs..see the documentaion')
      return(0)
    }
    h2 = c(0)
  }
  slices = 100
  if(!missing(numSlices)) {
    slices = numSlices 
  }
  n = 4000000
  if (!missing(N)) {
    n = N
  }
  
  construct_copula = TRUE
  
  # set flag for verbose mode
  verbose=FALSE
  
  #---------------------#
  round_only = FALSE
  rotate_only = FALSE
  W = 0
  e = 0
  internalpoint = c(0)
  Gaussian = FALSE
  win_len = 0
  NN = 0
  C = 0
  ratio = 0
  frac = 0
  ballwalk = FALSE
  delta =0
  vpoly = FALSE
  Zono = FALSE
  exact_zono = FALSE
  sample_only = FALSE
  var = 0
  coord = TRUE
  rounding = FALSE
  gen_only = FALSE
  Vpoly_gen = FALSE
  kind_gen = -1
  m_gen = 0
  dim_gen = 0
  exact_zono = FALSE
  ball_only = FALSE
  sam_simplex = FALSE
  sam_can_simplex = FALSE
  sam_arb_simplex = FALSE
  sam_ball = FALSE
  sam_sphere = FALSE
  sliceSimplex = FALSE
  #-------------------#
  
  Mat = vol_R(E, W, e, internalpoint, Gaussian, win_len, NN, C, ratio, frac,
                    ballwalk, delta, vpoly, Zono, exact_zono, gen_only, Vpoly_gen,
                    kind_gen, dim_gen, m_gen, round_only, rotate_only, ball_only,
                    sample_only, sam_simplex, sam_can_simplex, sam_arb_simplex, sam_ball,
                    sam_sphere, n, var, construct_copula, h1, h2, slices, sliceSimplex,
                    coord, rounding, verbose)
  
  return(Mat)
  
  
}