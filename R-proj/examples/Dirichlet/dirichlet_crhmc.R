# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2023 Vissarion Fisikopoulos
# Copyright (c) 2018-2023 Apostolos Chalkis
# Copyright (c) 2020-2023 Elias Tsigaridas

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for sampling from the Dirichlet distribution with CRHMC

# Import required libraries
library(volesti)
library(Matrix)

Aeq = matrix(c(1, 1, 1), ncol=3, nrow=1, byrow=TRUE)
Aeq = as( Aeq, 'dgCMatrix' )
beq = matrix(c(1), ncol=1, nrow=1, byrow=TRUE)
Aineq = matrix(, ncol=3, nrow=0, byrow=TRUE)
Aineq = as( Aineq, 'dgCMatrix' )
bineq = matrix(, ncol=0, nrow=0, byrow=TRUE)
lb = c(0,0,0)
ub = c(1, 1 ,1)

# Create domain of truncation
P <- volesti::sparse_constraint_problem$new(Aineq, bineq, Aeq, beq, lb, ub)

a_vec = c(2,3,4)
# Negative log-probability oracle
f <- function(x) (-((a_vec-1)%*%log(x))[1])

# Negative log-probability gradient oracle
grad_f <- function(x) (-(a_vec-1)/x)

n_samples <- 5000
n_burns <- n_samples / 2
pts <- sample_points(P, n = n_samples,
                     random_walk = list("walk" = "CRHMC", 
                     "nburns" = n_burns, "walk_length" = 1, 
                     "solver" = "implicit_midpoint"), 
                     distribution = list("density" = "logconcave", 
                     "negative_logprob" = f, 
                     "negative_logprob_gradient" = grad_f))