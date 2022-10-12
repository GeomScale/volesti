This is an R package based on [volesti](https://github.com/GeomScale/volume_approximation) library that serves as supplementary code to sample uniformly distributed portfolios

###  Dependencies

To run the code you need `R >= 3.6.3` and you have to install the following `R` packages:  

1. `volesti` dependencies (see the DESCRIPTION file in folder `root/R-proj`)  
2. `Rmosek` first complete the (a) [mosek installation guide](https://docs.mosek.com/9.2/install/installation.html) and then proceed to (b) [Rmosek installation guide](https://docs.mosek.com/9.2/rmosek/install-interface.html)  
3. `pracma`
4. `Matrix`
5. `R.matlab`

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
inner_ball = get_max_inner_ball(result_list$HP_rounded$A, result_list$HP_rounded$b)
more_samples = sample_points(result_list$HP_rounded, n = N, 
                             random_walk = list("walk" = "aBiW", "walk_length" = 1,
                                                "starting_point"=inner_ball$center, 
                                                "L" = 6*sqrt(result_list$HP_rounded$dimension)*inner_ball$radius))

## map the points back to P0
samples_in_P0 = result_list$T %*% more_samples + 
  kronecker(matrix(1, 1, N), matrix(result_list$T_shift, ncol = 1))

## compute the portfolios
more_random_portfolios = result_list$N %*% samples_in_P0 + 
  kronecker(matrix(1, 1, N), matrix(result_list$N_shift, ncol = 1))
```
