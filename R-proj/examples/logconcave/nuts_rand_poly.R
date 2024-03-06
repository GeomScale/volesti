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

dimension <- 50
facets <- 200

# Create domain of truncation
H <- gen_rand_hpoly(dimension, facets, seed = 15)

# Rounding
Tr <- rounding(H, seed = 127)

P <- Hpolytope$new(A = Tr$Mat[1:nrow(Tr$Mat), 2:ncol(Tr$Mat)], b = Tr$Mat[,1])

x_min = matrix(0, dimension, 1)

# Warm start point from truncated Gaussian
warm_start <- sample_points(P, n = 1, random_walk = list("nburns" = 5000), distribution = list("density" = "gaussian", "variance" = 1/2, "mode" = x_min))

# Sample points
n_samples <- 20000

samples <- sample_points(P, n = n_samples, random_walk = list("walk" = "NUTS", "solver" = "leapfrog", "starting_point" = warm_start[,1]),
                         distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f))

# Plot histogram
hist(samples[1,], probability=TRUE, breaks = 100)

psrfs <- psrf_univariate(samples)
n_ess <- ess(samples)

