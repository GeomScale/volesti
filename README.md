This is an R package based on [volesti](https://github.com/GeomScale/volume_approximation) library that serves as supplementary code to sample uniformly distributed portfolios

###  Dependencies

To run the code you need `R >= 3.6.3` and you have to install the following `R` packages:  

1. `volesti` dependencies (see the DESCRIPTION file in folder `root/R-proj`)  
2. `quadprog` 
3. `pracma`
4. `Matrix`

###  Installation

To install volesti of the current branch, in folder `/root/R-proj` run the following command in `R`:  
```r
Rcpp::compileAttributes()  
devtools::install()  
library(volesti)  
```

### Sample uniformly distributed portfolios  

Request effective sample size ess = 1000

```r
result_list = samples_uniform_portfolios(A=A, b=b, Aeq=Aeq, beq=beq, ess = 1000)
```

### Results  

result_list$HP_rounded : The full dimensional polytope of the last phase in MMCS method (use this to sample more portfolios).  
result_list$samples : The samples in the full dimensional space.  
result_list$random_portfolios : The uniformly distributed portfolios generated.  
result_list$N : The matrix of the null space.  
result_list$N_shift : Shift for the full dimensional polytope.  
result_list$T : The transformation matrix to map samples from the last phase to the initial phase.  
result_list$T_shift : The shift of the previous transformation.  
result_list$minWeights : The minimum value of each weight.  
result_list$maxWeights : The maximum value of each weight.  
result_list$run_time : The total run-time.  

To sample more portfolios, that is how you can use the last phase,  

```r
N = 5000 # sample 5000 more portfolios
more_random_portfolios = sample_from_last_polytope(result_list, N)
```
