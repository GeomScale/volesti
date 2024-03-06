# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou
# Copyright (c) 2022-2022 Ioannis Iakovidis

# Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of Code 2022 program.

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for using the logconcave sampling methods

# Import required libraries
library(volesti)

# Sampling from uniform density example

logconcave_sample<- function(P,distribution, n_samples ,n_burns){
  if (distribution == "uniform"){
    f <- function(x) (0)
    grad_f <- function(x) (0)
    L=1
    m=1
    pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "CRHMC", "nburns" = n_burns, "walk_length" = 1, "solver" = "implicit_midpoint"), distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f, "L_" = L, "m" = m))
    return(max(psrf_univariate(pts, "interval")))
  }
  else if(distribution == "gaussian"){
    pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "CRHMC", "nburns" = n_burns, "walk_length" = 1, "solver" = "implicit_midpoint"), distribution = list("density" = "logconcave", "variance"=8))
    return(max(psrf_univariate(pts, "interval")))
  }
}

for (i in 1:2) {

  if (i==1) {
    distribution = 'gaussian'
    cat("Gaussian ")
  } else {
    distribution = 'uniform'
    cat("Uniform ")
  }

  P = gen_simplex(10, 'H')
  psrf = logconcave_sample(P,distribution,5000,2000)
  cat("psrf = ")
  cat(psrf)
  cat("\n")
}
