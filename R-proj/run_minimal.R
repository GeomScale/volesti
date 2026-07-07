library(volesti)
n <- 5                                  # number of assets
set.seed(42)                            # reproducible random data
A <- matrix(rnorm(n * 3), 3, n)         # random factor loadings
sigma <- cov2cor(t(A) %*% A + diag(0.1, n))  # pos-def, scaled to correlations

c <- 0.3                                  # target portfolio variance

# Sample portfolios whose variance equals c = 0.3
result <- sample_ptfs_constant_volatility(
  sigma,                                # covariance matrix
  c    = c,                             # target portfolio variance
  M    = 2000                           # points per random walk
)

samples <- result$overall_samples[[1]]  # matrix: rows = assets, cols = portfolios
sum(abs(colSums(samples) - 1))          # long-only simplex constraint met
# check that all sampled portfolios have variance equal to c up to numerical precision
sum(abs(diag(t(samples) %*% sigma %*% samples) - c)) 