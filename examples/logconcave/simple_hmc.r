# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou

# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for using the logconcave sampling methods

# Import required libraries
library(ggplot2)
library(volesti)

# Sampling from logconcave density example

# Helper function
norm_vec <- function(x) sqrt(sum(x^2))

# Negative log-probability oracle
f <- function(x) (norm_vec(x)^2 + sum(x))

# Negative log-probability gradient oracle
grad_f <- function(x) (2 * x + 1)

# Interval [-1, 1]
A = matrix(c(1, -1), ncol=1, nrow=2, byrow=TRUE)
b = c(2,1)

# Create domain of truncation
P <- volesti::Hpolytope$new(A, b)

# Mode of logconcave density
x_min <- c(-0.5)

# Smoothness and strong-convexity
L <- 2
m <- 2

# Warm start point from truncated Gaussian
warm_start <- sample_points(P, n = 1, random_walk = list("nburns" = 5000), distribution = list("density" = "gaussian", "variance" = 1/L, "mode" = x_min))

# Sample points
n_samples <- 20000
n_burns <- n_samples / 2

pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "HMC", "step_size" = 0.3, "nburns" = n_burns, "walk_length" = 3, "solver" = "leapfrog", "starting_point" = warm_start[,1]), distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f, "L_" = L, "m" = m))

# Plot histogram
hist(pts, probability=TRUE, breaks = 100)

cat("Sample mean is: ")
sample_mean <- mean(pts)
cat(sample_mean)
cat("\n")
cat("Sample variance is: ")
sample_variance <- mean((pts - sample_mean)^2)
cat(sample_variance)
