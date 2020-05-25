#' A function to over-approximate a zonotope with PCA method and to evaluate the approximation by computing a ratio of fitness.
#' 
#' For the evaluation of the PCA method the exact volume of the approximation body is computed and the volume of the input zonotope is computed by CoolingBodies algorithm (BAN). The ratio of fitness is \eqn{R=}.
#' 
#' @param Z A zonotope.
#' @param fit_ratio Optional. A boolean parameter to request the computation of the ratio of fitness.
#' @param settings Optional. A list that declares the values of the parameters of CB algorithm as follows:
#' \itemize{
#' \item{\code{error} }{ A numeric value to set the upper bound for the approximation error. The default value is \eqn{1} for \code{'SOB'} and \eqn{0.1} otherwise.}
#' \item{\code{walk_length} }{ An integer to set the number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for \code{'SOB'} and \eqn{1} otherwise.}
#' \item{\code{win_len} }{ The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
#' \item{\code{hpoly} }{ A boolean parameter to use H-polytopes in MMC of CB algorithm. The default value is \code{FALSE}.}
#' }
#' @param seed Optional. A fixed seed for the number generator.
#' 
#' @return A list that contains the approximation body in H-representation and the ratio of fitness
#' 
#' @examples 
#' # over-approximate a 2-dimensional zonotope with 10 generators and compute the ratio of fitness
#' Z = gen_rand_zonotope(2,10)
#' retList = zonotope_approximation(Z = Z, fit_ratio = TRUE)
#' 
#' @export
zonotope_approximation <- function(Z, fit_ratio = NULL, settings = NULL, seed = NULL){
  
  ret_list = zono_approx(Z, fit_ratio, settings, seed)
  
  Mat = ret_list$Mat
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  PP = list("P" = Hpolytope$new(A,b), "fit_ratio" = ret_list$fit_ratio)
  
  return(PP)
  
}

