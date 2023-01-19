# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou
# Copyright (c) 2022-2022 Ioannis Iakovidis

# Contributed and/or modified by Ioannis Iakovidis, as part of Google Summer of Code 2022 program.

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

# Sample points
n_samples <- 80000
n_burns <- n_samples / 2
cat("---Sampling without hessian\n")
pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "CRHMC", "step_size" = 0.3, "nburns" = n_burns, "walk_length" = 1, "solver" = "implicit_midpoint"), distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f, "L_" = L, "m" = m))
jpeg("histogram_without_hessian.jpg")
# Plot histogram
hist(pts, probability=TRUE, breaks = 100)

cat("Sample mean is: ")
sample_mean <- mean(pts)
cat(sample_mean)
cat("\n")
cat("Sample variance is: ")
sample_variance <- mean((pts - sample_mean)^2)
cat(sample_variance)
cat("\n")
invisible(capture.output(dev.off()))

# Negative log-probability hessian oracle
hess_f <- function(x) (2)
cat("---Sampling with hessian\n")
pts <- sample_points(P, n = n_samples, random_walk = list("walk" = "CRHMC", "step_size" = 0.3, "nburns" = n_burns, "walk_length" = 1, "solver" = "implicit_midpoint"), distribution = list("density" = "logconcave", "negative_logprob" = f, "negative_logprob_gradient" = grad_f,"negative_logprob_hessian" = hess_f, "L_" = L, "m" = m))
jpeg("histogram_with_hessian.jpg")
# Plot histogram
hist(pts, probability=TRUE, breaks = 100)

cat("Sample mean is: ")
sample_mean <- mean(pts)
cat(sample_mean)
cat("\n")
cat("Sample variance is: ")
sample_variance <- mean((pts - sample_mean)^2)
cat(sample_variance)
cat("\n")
invisible(capture.output(dev.off()))

walk="CRHMC"
library(Matrix)
bineq=matrix(c(10,10,10,10,10), nrow=5, ncol=1, byrow=TRUE)
Aineq = matrix(c(1,0,-0.25,-1,2.5,1,0.4,-1,-0.9,0.5), nrow=5, ncol=2, byrow = TRUE)
Aineq = as( Aineq, 'dgCMatrix' )
beq=matrix(,nrow=0, ncol=1, byrow=TRUE)
Aeq = matrix(, nrow=0, ncol=2, byrow = TRUE)
Aeq=as( Aeq, 'dgCMatrix' )
lb=-100000*c(1,1);
ub=100000*c(1,1);
cat("---Sampling the normal distribution in a pentagon\n")
P <- volesti::sparse_constraint_problem$new(Aineq, bineq,Aeq, beq, lb, ub)
points <- sample_points(P, n = n_samples, random_walk = list("walk" = "CRHMC", "step_size" = 0.3, "nburns" = n_burns, "walk_length" = 1, "solver" = "implicit_midpoint"), distribution = list("density" = "logconcave", "variance" = 8))
jpeg("pentagon.jpg")
plot(ggplot(data.frame( x=points[1,], y=points[2,] )) +
geom_point( aes(x=x, y=y, color=walk)) + coord_fixed(xlim = c(-15,15),
ylim = c(-15,15)) + ggtitle(sprintf("Sampling a random pentagon with walk %s", walk)))
invisible(capture.output(dev.off()))
write.table(points, file="pentagon.txt", row.names=FALSE, col.names=FALSE)
