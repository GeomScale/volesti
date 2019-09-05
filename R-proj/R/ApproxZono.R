#' A function to over-approximate a zonotope with PCA method and to evaluate the approximation by computing a ratio of fitness.
#' 
#' For the evaluation of the PCA method the exact volume of the approximation body is computed and the volume of the input zonotope is computed by CoolingBodies algorithm (BAN). The ratio of fitness is \eqn{R=}.
#' 
#' @param Z A zonotope.
#' @param fit_ratio Optional. A boolean parameter to request the computation of the ratio of fitness.
#' @param walk_length Optional. The number of the steps for the random walk. The default value is \eqn{\lfloor 10 + d/10\rfloor} for SequenceOfBalls and \eqn{1} for CoolingGaussian.
#' @param error Optional. Declare the upper bound for the approximation error. The default value is \eqn{1} for SequenceOfBalls and \eqn{0.1} for CoolingGaussian.
#' @param inner_ball Optional. A \eqn{d+1} vector that contains an inner ball. The first \eqn{d} coordinates corresponds to the center and the last one to the radius of the ball. If it is not given then for H-polytopes the Chebychev ball is computed, for V-polytopes \eqn{d+1} vertices are picked randomly and the Chebychev ball of the defined simplex is computed. For a zonotope that is defined by the Minkowski sum of \eqn{m} segments we compute the maximal \eqn{r} s.t.: \eqn{re_i\in Z} for all \eqn{i=1,\dots ,d}, then the ball centered at the origin with radius \eqn{r/\sqrt{d}} is an inscribed ball.
#' @param random_walk Optional. A string that declares the random walk method: a) \code{'CDHR'} for Coordinate Directions Hit-and-Run, b) \code{'RDHR'} for Random Directions Hit-and-Run or c) \code{'BW'} for Ball Walk. The default walk is \code{'CDHR'}.
#' @param rounding Optional. A boolean parameter for rounding. The default value is \code{FALSE}.
#' @param parameters Optional. A list for the parameters of the algorithms:
#' \itemize{
#'  \item{\code{Window} }{ The length of the sliding window for CG algorithm. The default value is \eqn{500+4dimension^2}.}
#'  \item{\code{BW_rad} }{ The radius for the ball walk. The default value is \eqn{4r/dimension}, where \eqn{r} is the radius of the inscribed ball of the polytope.}
#'  \item{\code{ub} }{ The lower bound for the ratios in MMC in BAN algorithm. The default value is \eqn{0.1}.}
#'  \item{\code{lb} }{ The upper bound for the ratios in MMC in BAN algorithm. The default value is \eqn{0.15}.}
#'  \item{\code{N} }{ An integer that controls the number of points \eqn{\nu N} generated in each convex body in annealing schedule.}
#'  \item{\code{nu} }{ The degrees of freedom for the t-student distribution in t-tests in BAN algorithm. The default value is \eqn{10}.}
#'  \item{\code{alpha} }{ The significance level for the t-tests in BAN algorithm. The default values is 0.2.}
#'  \item{\code{prob} }{ The probability is used for the empirical confidence interval in ratio estimation of BAN algorithm. The default value is \eqn{0.75}.}
#'  \item{\code{hpoly} }{ A boolean parameter to use H-polytopes in MMC of BAN algorithm. The default value is \code{FALSE}.}
#'  \item{\code{minmaxW} }{ A boolean parameter to use the sliding window with a minmax values stopping criterion.}
#' }
#' 
#' @return A list that contains the approximation body in H-representation and the ratio of fitness
#' 
#' @examples 
#' # over-approximate a 2-dimensional zonotope with 10 generators and compute the ratio of fitness
#' Z = GenZonotope(2,10)
#' retList = ApproxZono(Z = Z, fit_ratio = TRUE)
#' 
#' @export
ApproxZono <- function(Z, fit_ratio = NULL, walk_length = NULL, error = NULL,
                        inner_ball = NULL, random_walk = NULL, rounding = NULL,
                        parameters = NULL){
  
  ret_list = zono_approx(Z, fit_ratio, walk_length, error, inner_ball, random_walk, rounding, parameters)
  
  Mat = ret_list$Mat
  
  # first column is the vector b
  b = Mat[,1]
  
  # remove first column
  A = Mat[,-c(1)]
  PP = list("P" = Hpolytope$new(A,b), "fit_ratio" = ret_list$fit_ratio)
  
  return(PP)
  
}