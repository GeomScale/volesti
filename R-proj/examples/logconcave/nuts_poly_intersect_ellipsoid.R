# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-223 Elias Tsigaridas

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

# Create domain of truncation
A = matrix(c(1,0,0,1,-1,0,0,-1), nrow=4, ncol=2, byrow=TRUE)
b = rep(1,4)
E = matrix(c(0.25, 0.75, 0.75, 3.25), nrow=2, ncol=2, byrow=TRUE)
EP = PolytopeIntersectEllipsoid$new(A, b, E)
dim = 2

x0 = matrix(0, dim, 1)

# Sample points
n_samples <- 2000

samples <- sample_points(EP, n = n_samples, random_walk = list("walk" = "NUTS", "solver" = "leapfrog", "starting_point" = x0),
                         distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f))

# Plot histogram
hist(samples[1,], probability=TRUE, breaks = 100)

psrfs <- psrf_univariate(samples)
n_ess <- ess(samples)


