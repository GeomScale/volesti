# Code for "Reassessing the Low-Volatility Anomaly through the Geometry of Portfolio Choice"

Companion code for the paper. Samples long-only constant-volatility portfolios by
geometric random walks on the intersection of the simplex with a volatility ellipsoid.

## Install

Requires the `volesti` R package. From the repo root:

```r
Rcpp::compileAttributes()
install.packages("R-proj", repos = NULL, type = "source")
```

Dependencies: `nloptr`, `MASS`, `Matrix`, `pracma`, `BiocManager`, `graph`, `Rcpp`, `RcppEigen`, `BH`.  
For `nloptr` on Ubuntu: `sudo apt install libnlopt-dev`.  
For `graph`: `BiocManager::install("graph")`.

## Quick start

```r
library(volesti)                        # load the package

# Build a synthetic 5-asset covariance matrix
n <- 5                                  # number of assets
set.seed(42)                            # reproducible random data
A <- matrix(rnorm(n * 3), 3, n)         # random factor loadings
sigma <- cov2cor(t(A) %*% A + diag(0.1, n))  # pos-def, scaled to correlations

# Sample portfolios whose variance equals c = 0.3
result <- sample_ptfs_constant_volatility(
  sigma,                                # covariance matrix
  c    = 0.3,                           # target portfolio variance
  M    = 2000                           # points per random walk
)

samples <- result$overall_samples[[1]]  # matrix: rows = assets, 
                                        # cols = portfolios
colSums(samples)                        # long-only simplex constraint
mean(diag(t(samples) %*% sigma %*% samples))  # ~0.3; target volatility
```

See `run_minimal_example.R` for a complete worked example.

## Reproducing the paper

The full pipeline requires covariance matrices from real market data
(`usa_covariance_matrices_small.rds`, `usa_covariance_matrices_large.rds`).
These are not bundled here - contact the authors. Once available, run:

```r
# Small universe (fewer assets)
source("run_real_data_small.R")

# Large universe (more assets) 
source("run_real_data_large.R")
```

Both scripts iterate over time windows and volatility levels, calling
`sample_ptfs_constant_volatility()` for each (sigma, c) pair, then save
results as `.rds` files.

## Core function

`sample_ptfs_constant_volatility(sigma, c, M, ignore_smallest_components = TRUE)`

- `sigma` — nxn covariance matrix
- `c` — target portfolio variance
- `M` — number of points to sample per walk
- Returns a list with `$overall_samples` (matrix, columns = portfolios) and `$num_verts`
