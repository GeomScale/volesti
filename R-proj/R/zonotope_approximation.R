#' A function to over-approximate a zonotope with PCA method and to evaluate the approximation by computing a ratio of fitness.
#' 
#' For the evaluation of the PCA method the exact volume of the approximation body is computed and the volume of the input zonotope is computed by CoolingBodies algorithm (BAN). The ratio of fitness is \eqn{R=}.
#' 
#' @param Z A zonotope.
#' @param fit_ratio Optional. A boolean parameter to request the computation of the ratio of fitness.
#' @param algo_parameters Optional. A list that declares the values of the parameters of CB algorithm as follows:
#' \itemize{
#' \item{\code{error} }{ A numeric value to set the upper bound for the approximation error. The default value is \eqn{1} for \code{'SOB'} and \eqn{0.1} otherwise.}
#' \item{\code{random_walk} }{ A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run, c) \code{'BaW'} for Ball Walk, or \code{'BiW'} for Billiard walk. The default walk is \code{'CDHR'} for H-polytopes and \code{'BiW'} for the other representations.}
#' \item{\code{walk_length} }{ An integer to set the number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for \code{'SOB'} and \eqn{1} otherwise.}
#' \item{\code{rounding} }{ A boolean parameter for rounding. The default value is \code{FALSE}.}
#' \item{\code{inner_ball} }{  A \eqn{d+1} numeric vector that contains an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an inscribed ball.}
#' \item{\code{len_win} }{ The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
#' \item{\code{BW_rad} }{ The radius for the ball walk. The default value is \eqn{4r/dimension}, where \eqn{r} is the radius of the inscribed ball of the polytope.}
#' \item{\code{ub} }{ The lower bound for the ratios in MMC in CB algorithm. The default value is \eqn{0.1}.}
#' \item{\code{lb} }{ The upper bound for the ratios in MMC in CB algorithm. The default value is \eqn{0.15}.}
#' \item{\code{N} }{ An integer that controls the number of points \eqn{\nu N} generated in each convex body in annealing schedule of algorithm CB.}
#' \item{\code{nu} }{ The degrees of freedom for the t-student distribution in t-tests in CB algorithm. The default value is \eqn{10}.}
#' \item{\code{alpha} }{ The significance level for the t-tests in CB algorithm. The default values is 0.2.}
#' \item{\code{prob} }{ The probability is used for the empirical confidence interval in ratio estimation of CB algorithm. The default value is \eqn{0.75}.}
#' \item{\code{hpoly} }{ A boolean parameter to use H-polytopes in MMC of CB algorithm. The default value is \code{FALSE}.}
#' \item{\code{minmaxW} }{ A boolean parameter to use the sliding window with a minmax values as a stopping criterion.}
#' \item{\code{L} }{The maximum length of the billiard trajectory.}
#' }
#' 
#' @return A list that contains the approximation body in H-representation and the ratio of fitness
#' 
#' @examples 
#' # over-approximate a 2-dimensional zonotope with 10 generators and compute the ratio of fitness
#' Z = gen_rand_zonotope(2,10)
#' retList = zonotope_approximation(Z = Z, fit_ratio = TRUE)
#' 
#' @export
zonotope_approximation <- function(Z, fit_ratio = NULL, algo_parameters = NULL){
  
  ret_list = zono_approx(Z, fit_ratio, algo_parameters)
  
  Mat = ret_list$Mat
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  PP = list("P" = Hpolytope$new(A,b), "fit_ratio" = ret_list$fit_ratio)
  
  return(PP)
  
}

