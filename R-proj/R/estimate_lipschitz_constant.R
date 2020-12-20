library(fitdistrplus)

dinvweibull_with_loc <- function (x, k, lambda, theta) {
  return ((k / lambda) * ((x - theta) / lambda)^(k - 1) * exp(- ((x - theta) / lambda)^k) * as.double(x >= 0))
}

pinvweibull_with_loc <- function (q, k, lambda, theta) {
  return ((1 - exp(-((q - theta) / lambda)^k)) * as.double(q >= 0))
}


#' Estimate the Lipschitz Constant of a function f
#'
#' @export
estimtate_lipschitz_constant <- function (f, P, n) {
  points = sample_points(P, n = 1000, random_walk = list("walk" = "BaW", "walk_length" = 1))  
  l = matrix(0, 1, n-1)
  
  for (i in seq(n-1)) {
    l[1, i] <- norm(as.matrix(f(points[,i+1])) - as.matrix(f(points[,i])), type='F') / norm(as.matrix(points[,i+1]) - as.matrix(points[,i]), type='F')
  }
  
  na_cols <- is.nan(l[1,])
  
  l <- as.matrix(l[1, !na_cols])
  
  return(max(l))
}