#' Inverse weibull distribution PDF with location parameter
#' @param x The argument of the PDF
#' @param k The shape parameter
#' @param lambda The scale parameter
#' @param theta The location parameter
#' 
#' @return The value of the PDF of an Inverse Weibull distribution with parameters k, lambda, theta evaluated at x
#' @export
dinvweibull_with_loc <- function (x, k, lambda, theta) {
  return ((k / lambda) * ((x - theta) / lambda)^(k - 1) * exp(- ((x - theta) / lambda)^k) * as.double(x >= 0))
}

#' Inverse weibull distribution CDF with location parameter
#' @param q The argument of the CDF
#' @param k The shape parameter
#' @param lambda The scale parameter
#' @param theta The location parameter
#' 
#' @return The value of the CDF of an Inverse Weibull distribution with parameters k, lambda, theta evaluated at q
#' @export
pinvweibull_with_loc <- function (q, k, lambda, theta) {
  return ((1 - exp(-((q - theta) / lambda)^k)) * as.double(q >= 0))
}


#' Estimate the Lipschitz Constant of a function f
#' 
#' @param f Function whose Lipschitz constant is to be estimated
#' @param P Domain of f (a convex polytope)
#' @param n Number of samples to take
#' 
#' The procedure draws n uniform samples from P and evaluates the Lipschitz
#' constant at subsequent samples (where the sampler moves to a new point),
#' It then returns the maximum observation
#' 
#' @return An estimate of the Lipschitz constant
#' 
#' @export 
estimtate_lipschitz_constant <- function (f, P, n) {
  points = volesti::sample_points(P, n = 1000, random_walk = list("walk" = "BaW", "walk_length" = 1))
  l = matrix(0, 1, n-1)

  for (i in seq(n-1)) {
    l[1, i] <- norm(as.matrix(f(points[,i+1])) - as.matrix(f(points[,i])), type='F') / norm(as.matrix(points[,i+1]) - as.matrix(points[,i]), type='F')
  }

  na_cols <- is.nan(l[1,])

  l <- as.matrix(l[1, !na_cols])
  
  # TODO Implement weibull distribution fitting method to estimate the Lipschitz constant
  return(max(l))
}