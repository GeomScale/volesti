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
#' @export
copula1 <-function(h1, h2, numSlices, N) {
  
  if(missing(h1) || missing(h2)) {
    stop('Only one family is given.')
  }
  
  if (length(h1)<2 || length(h2)<2) {
    stop("The dimension of the two families has to be greater than 2.")
  }
  
  if (length(h1)!=length(h2)){
    stop('The two families are defined for different dimensions.')
  }
  
  slices = 100
  if(!missing(numSlices)) {
    slices = numSlices
    if(slices<=0) {
      stop("Number of slices has to be a positive integer.")
    }
  }
  n = 4000000
  if (!missing(N)) {
    n = N
    if(n<=0) {
      stop("Number of sampled points has to be a positive integer.")
    }
  }
  
  Mat = copula_hyps(h1, h2, slices, n)
  return(Mat)
}