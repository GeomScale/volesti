pinv <- function(A, eps=1e-8){
  L <- svd(A)
  d <- L$d
  i <- abs(d) > eps
  d[i] <- 1/d[i]
  L$v %*% diag(d, nrow=length(d)) %*% t(L$u)
}