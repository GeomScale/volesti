# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou

# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for sampling from a Generalized Hyperbolic density

# Import required libraries
library(ggplot2)
library(volesti)
library(numDeriv)
library(GeneralizedHyperbolic)

A = matrix(c(1, -1), ncol=1, nrow=2, byrow=TRUE)
b = c(4,4)

f <- function(x) (-log(dghyp(x)))
grad_f <- function(x) (-ddghyp(x)/dghyp(x))

x_min = matrix(0, 1, 1)

# Create domain of truncation
P <- volesti::Hpolytope$new(A, b)

# Smoothness and strong-convexity
L <- estimtate_lipschitz_constant(grad_f, P, 1000)
m <- L

# Warm start point from truncated Gaussian
warm_start <- sample_points(P, n = 1, random_walk = list("nburns" = 5000), distribution = list("density" = "gaussian", "variance" = 1/L, "mode" = x_min))

# Sample points
n_samples <- 10000
n_burns <- n_samples / 2

pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "HMC", "step_size" = 0.5, "nburns" = n_burns, "walk_length" = 1, "solver" = "leapfrog", "starting_point" = warm_start[,1]), distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f, "L_" = L, "m" = m))

# Plot histogram
hist(pts, 
     probability=TRUE, 
     breaks = 100, 
     border="blue",
     main="Genrealized Hyperbolic Density with lambda = 1, alpha = 1, beta = 0, delta = 1, mu = 0",
     xlab="Samples",
     ylab="Density"
)

cat("Sample mean is: ")
sample_mean <- mean(pts)
cat(sample_mean)
cat("\n")
cat("Sample variance is: ")
sample_variance <- mean((pts - sample_mean)^2)
cat(sample_variance)

n_ess = min(ess(pts))
psrf = max(psrf_univariate(pts))

cat("\nEffective sample size: ", n_ess, append=TRUE)
cat("\nPSRF: ", psrf, append=TRUE)
