#' Compute an indicator for each time period that describes the state of a market.
#'
#' Given a matrix that contains row-wise the assets' returns and a sliding window \code{win_length}, this function computes an approximation of the joint distribution (copula, e.g. see \url{https://en.wikipedia.org/wiki/Copula_(probability_theory)}) between portfolios' return and volatility in each time period defined by \code{win_len}. 
#' For each copula it computes an indicator: If the indicator is large it corresponds to a crisis period and if it is small it corresponds to a normal period. 
#' In particular, the periods over which the indicator is greater than 1 for more than 60 consecutive sliding windows are warnings and for more than 100 are crisis. The sliding window is shifted by one day.
#'
#' @param returns A \eqn{d}-dimensional vector that describes the direction of the first family of parallel hyperplanes.
#' @param parameters A list to set a parameterization.
#' \itemize{
#' \item{win_length }{ The length of the sliding window. The default value is 60.}
#' \item{m } { The number of slices for the copula. The default value is 100.}
#' \item{n }{ The number of points to sample. The default value is \eqn{5\cdot 10^5}.}
#' \item{nwarning }{ The number of consecutive indicators larger than 1 required to declare a warning period. The default value is 60.}
#' \item{ncrisis }{ The number of consecutive indicators larger than 1 required to declare a crisis period. The default value is 100.}
#' \item{seed }{ A fixed seed for the number generator.}
#' }
#'
#' @references \cite{L. Cales, A. Chalkis, I.Z. Emiris, V. Fisikopoulos,
#' \dQuote{Practical volume computation of structured convex bodies, and an application to modeling portfolio dependencies and financial crises,} \emph{Proc. of Symposium on Computational Geometry, Budapest, Hungary,} 2018.}
#'
#' @return A list that contains the indicators and the corresponding vector that label each time period with respect to the market state: a) normal, b) crisis, c) warning.
#'
#' @examples 
#' # simple example on random asset returns
#' asset_returns = replicate(10, rnorm(14))
#' market_states_and_indicators = compute_indicators(asset_returns,
#'     parameters = list("win_length" = 10, "m" = 10, "n" = 10000, "nwarning" = 2, "ncrisis" = 3))
#'
#' @export
#' @useDynLib volesti, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp loadModule
#' @importFrom "utils" "read.csv"
#' @importFrom "stats" "cov"
#' @importFrom "methods" "new"
compute_indicators <- function(returns, parameters = list("win_length" = 60, "m" = 100, "n" = 500000, "nwarning" = 60, "ncrisis" = 100)) {
  
  win_length = 60
  if (!is.null(parameters$win_length)) {
    win_length = parameters$win_length
  }
  m=100
  if (!is.null(parameters$m)){
    m = parameters$m
  }
  n = 500000
  if (!is.null(parameters$n)){
    n = parameters$n
  }
  nwarning = 60
  if (!is.null(parameters$nwarning)) {
    nwarning = parameters$nwarning
  }
  ncrisis = 100
  if (!is.null(parameters$ncrisis)) {
    ncrisis = parameters$ncrisis
  }
  seed = NULL
  if (!is.null(parameters$seed)) {
    seed = parameters$seed
  }
  
  nrows = dim(returns)[1]
  nassets = dim(returns)[2]
  wl = win_length-1
  
  indicators = c()
  for (i in 1:(nrows-wl)) {
    
    Win=i:(i+wl)
    E = cov(returns[Win,])
    
    compRet = rep(1,nassets)
    for (j in 1:nassets) {
      for (k in Win) {
        compRet[j] = compRet[j] * (1 + returns[k, j])
      }
      compRet[j] = compRet[j] - 1
    }
    
    cop = copula(r1 = compRet, sigma = E, m = m, n = n, seed = seed)
    blue_mass = 0
    red_mass = 0
    
    for (row in 1:m) {
      for (col in 1:m) {
        if (row-col<=0.2*m && row-col>=-0.2*m) {
          if (row+col<0.8*m || row+col>1.2*m) {
            red_mass = red_mass + cop[row,col]
          }
        } else {
          if (row+col>=0.8*m+1 && row+col<=1.2*m+1) {
            blue_mass = blue_mass + cop[row,col]
          }
        }
      }
    }
    indicators = c(indicators, blue_mass / red_mass)
  }

  N = length(indicators)

  index = 0
  set_index = FALSE
  col = rep("normal", N)

  for (i in 1:N) {
  
    if(indicators[i]>1 && !set_index){
      index = i
      set_index = TRUE
    } else if (indicators[i]<1) {
      if(set_index){
        if(i-index > nwarning-1 && i-index <= ncrisis-1){
          col[index:(i-1)] = "warning"
        } else if(i-index > ncrisis-1) {
          col[index:(i-1)] = "crisis"
        }
      }
      set_index = FALSE
    }
  }
  if(set_index){
    if(N-index+1 > nwarning-1 && N-index+1 <= ncrisis-1){
      col[index:i] = "warning"
    } else if(N-index+1 > ncrisis-1) {
      col[index:i] = "crisis"
    }
  }
  
  return(list("indicators" = indicators, market_states = col))
  
}
