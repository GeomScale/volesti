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
library(R.matlab)

# Sampling from logconcave density example

# Helper function
norm_vec <- function(x) sqrt(sum(x^2))


# Load polytopes from mat file
metabolic_polytope_mat <- readMat('./data/polytope_e_coli.mat')
A <- as.matrix(metabolic_polytope_mat$polytope[[1]])
b <- as.matrix(metabolic_polytope_mat$polytope[[2]])
center <- as.matrix(metabolic_polytope_mat$polytope[[3]])
radius <- as.numeric(metabolic_polytope_mat$polytope[[4]])
sigma <- 1
dimension <- dim(A)[2]


# Negative log-probability oracle
f <- function(x) (norm_vec(x)^2 / (2 * sigma^2))

# Negative log-probability gradient oracle
grad_f <- function(x) (x / sigma^2)


# Smoothness and strong-convexity
L <- 1 / sigma^2
m <- 1 / sigma^2

# Center polytope
b_new <- b - A %*% center

# Create volesti polytope
P <- Hpolytope$new(A = A, b = c(b_new))

# Rounding
#Tr <- rounding(H)

#P <- Hpolytope$new(A = Tr$Mat[1:nrow(Tr$Mat), 2:ncol(Tr$Mat)], b = Tr$Mat[,1])

# Center is origin (after shift)
x_min = matrix(0, dimension, 1)

# Generate samples with HNR
start_time <- Sys.time()
rdhr_samples <- sample_points(P, n = 10, random_walk = list("walk" = "RDHR", "nburns" = 10, "walk_length" = 1), distribution = list("density" = "gaussian", "variance" = 1/L, "mode" = x_min))
end_time <- Sys.time()

# Calculate Effective Sample size
rdhr_ess = ess(rdhr_samples)
min_ess <- min(rdhr_ess)

# Calculate PSRF
rdhr_psrfs = psrf_univariate(rdhr_samples)
max_psrf = max(rdhr_psrfs)
elapsed_time <- end_time - start_time

# Print results
cat('Min Effective Sample Size: ')
cat(min_ess)
cat('\n')
cat('Maximum PSRF: ')
cat(max_psrf)
cat('\n')
cat('Time per independent sample: ')
cat(elapsed_time / min_ess)
cat('sec')

outfile <- '/home/marios/samples_hnr_iAB_PLT_283.txt'

write.table(rdhr_samples, file=outfile, row.names=FALSE, col.names=FALSE)

start_time <- Sys.time()
hmc_samples <- sample_points(P, n = 10, random_walk = list("walk" = "HMC", "step_size" = 0.07, "nburns" = 10, "walk_length" = 30, "solver" = "leapfrog", "starting_point" = rdhr_samples[, ncol(rdhr_samples)]), distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f, "L_" = L, "m" = m))
end_time <- Sys.time()

# Calculate Effective Sample size
hmc_ess = ess(hmc_samples)
min_ess <- min(hmc_ess)

# Calculate PSRF
hmc_psrfs = psrf_univariate(hmc_samples)
max_psrf = max(hmc_psrfs)
elapsed_time <- end_time - start_time

# Print results
cat('HMC\n')
cat('Min Effective Sample Size: ')
cat(min_ess)
cat('\n')
cat('Maximum PSRF: ')
cat(max_psrf)
cat('\n')
cat('Time per independent sample: ')
cat(elapsed_time / min_ess)
cat('sec')
