# --------------------------------------------------------------------------
# Quadratic Inverse Shrinkage (NL) Ledoit and Wolf (2021)
# --------------------------------------------------------------------------
cov.qis <- function(Y, k = -1) 
{
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  sample <- (t(Y) %*% Y) / n    # sample covariance matrix    
  spectral <- eigen(sample)    # spectral decompositon
  lambda <- spectral$values[p:1]    # sort eigenvalues in ascending order
  u <- spectral$vectors[,p:1]    # eigenvectors follow their eigenvalues
  h <- min(c^2, 1/c^2)^0.35 / p^0.35    # smoothing parameter
  invlambda <- 1 / lambda[max(1, p-n+1):p]    # inverse of non-null eigenvalues   
  #Lj <- rep.row(invlambda, min(p, n))    # like 1 / lambda_j
  Lj <- matrix(rep(invlambda,each=min(p, n)),nrow=min(p, n))    # like 1 / lambda_j
  Lj.i <- Lj - t(Lj)    # like (1 / lambda_j) - (1 / lambda_i)
  theta <- rowMeans(Lj * Lj.i / (Lj.i^2 + h^2 * Lj^2))    # smoothed Stein shrinker
  Htheta <- rowMeans(Lj * (h * Lj) / (Lj.i^2 + h^2 * Lj^2)) # its conjugate
  Atheta2 <- theta^2 + Htheta^2    # its squared amplitude
  if (p <= n)    # case where sample covariance matrix is not singular
    delta <- 1 / ((1 - c)^2 * invlambda + 2 * c * (1 - c) * invlambda * theta +
                    c^2 * invlambda * Atheta2)           # optimally shrunk eigenvalues
  else {    # case where sample covariance matrix is singular
    delta0 <- 1 / ((c - 1) * mean(invlambda))     # shrinkage of null eigenvalues
    delta <- c(rep(delta0, p - n), 1 / (invlambda * Atheta2));
  }
  deltaQIS <- delta * (sum(lambda) / sum(delta))    # preserve trace
  sigmahat <- u %*% diag(deltaQIS) %*% t(u)    #reconstruct covariance matrix
  
  return( sigmahat )
}