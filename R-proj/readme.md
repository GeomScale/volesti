## Dependencies

R packages: `nloptr`, `MASS`, `Matrix`, `pracma`, `ggm`, `BiocManager`, `graph`.  

To install `graph` use `BiocManager::install("graph")`.  

## Installation

Run the followings:  

`Rcpp::compileAttributes()`  
`R CMD INSTALL --no-multiarch --with-keep.source R-proj`  

Or use the command `Build -> Install and Restart` in Rstudio.  

## How to use the code

To run the pipeline using data of assets returns, see the script `run_data.R`.  

To sample `M` long-only portfolios with volatily `c`, when the covariance matrix is `sigma`, run:  

```
Ptfs = sample_ptfs_constant_volatility(sigma, c, M)
```

