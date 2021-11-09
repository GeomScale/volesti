## Dependencies

R packages: `MASS`, `Matrix`, `pracma`, `ggm`, `BiocManager`, `graph`.  

To install `graph` use `BiocManager::install("graph")`.  

## Installation

Run the followings:  
  
`Rcpp::compileAttributes()`  
`R CMD INSTALL --no-multiarch --with-keep.source R-proj`  

## Use the code

See the script `run_data.R`.  
  
To sample `M` long-only portfolios with volatily `c`, when the covariance matrix is `sigma`, run:  
  
```
parameters = get_parameters(n-1)  # n is the number of asets
Ptfs = sample_ptfs_constant_volatility(sigma, c, M, parameters)
```


