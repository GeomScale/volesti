# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2023 Vissarion Fisikopoulos
# Copyright (c) 2018-2023 Apostolos Chalkis
# Copyright (c) 2020-2023 Elias Tsigaridas

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for sampling from the Dirichlet distribution with CRHMC

# Import required libraries
library(volesti)
library(Matrix)

n = 3
Aeq = matrix(rep(1,n), ncol=n, nrow=1, byrow=TRUE)
beq = c(1)
A = -diag(n)
b = rep(0, n)

# ellipsoid x'Ex <= 1
E = matrix(c(7.8791, 3.6018, -6.9191, 3.6018, 2.2995,
             -4.4496, -6.9191, -4.4496, 11.7993), 
           nrow = 3, ncol = 3, byrow = TRUE)

x0 = c(0.2455, 0.4287, 0.3258) #interior point

# Projecting everything onto the null space of Aeq
N_trans = pracma::nullspace(Aeq)
N_shift = x0

b = b - A %*% N_shift
A = A %*% N_trans
p0 = rep(0, n-1) #interior point of the projected point

#the new matrix of the projected ellipsoid
Et = t(N_trans) %*% E %*% N_trans
#the new center of the projected ellipsoid
center = -MASS::ginv(Et) %*% (t(N_trans) %*% E) %*% N_shift
#the new offset of the projeced ellipsoid
R = 1 + t(center)%*%Et%*%center - t(N_shift)%*%E%*%N_shift
R = R[1]
#normalizing the offset
Et = Et / R

#we shift by the ellipsoid's center so that it is given by x'Etx <= 1 
b = b - A %*% center # shift ellipsoid's center to origin
# interior point in the projected body
p0 = p0 - center

#now, from a point y, in the projected body, we get
#the point oin the initial body as follows: x = N_trans%*%(y+center) + x0

# Create domain of truncation
P <- volesti::PolytopeIntersectEllipsoid$new(A, b, Et)

a_vec = c(2,3,4)
# Transformed negative log-probability oracle
f <- function(y) (-((a_vec-1)%*%log(N_trans%*%(y+center) + x0))[1])

# Transformed negative log-probability gradient oracle
grad_f <- function(y) (-(a_vec-1)/(N_trans%*%(y+center) + x0))

n_samples <- 5000
n_burns <- n_samples / 2
samples <- volesti::sample_points(P, n = n_samples, 
                                  random_walk = list("walk" = "NUTS", "solver" = "leapfrog", 
                                                     "starting_point" = p0),
                                  distribution = list("density" = "logconcave", "negative_logprob" = f, 
                                                      "negative_logprob_gradient" = grad_f))

M = dim(samples)[2]
# samples in the initial body following Dirichlet distribution
random_portfolios = N_trans %*% (samples + kronecker(matrix(1, 1, M), matrix(center, ncol = 1))) + 
  kronecker(matrix(1, 1, M), matrix(N_shift, ncol = 1))

######################
## UNIFORM SAMPLING ##
######################
# it's better to use Billiard walk instead of a dirichlet with ones

samples <- volesti::sample_points(P, n = n_samples, 
                                  random_walk = list("walk" = "aBiW", "starting_point" = p0))

M = dim(samples)[2]
# samples in the initial body following uniform distribution
random_portfolios = N_trans %*% (samples + kronecker(matrix(1, 1, M), matrix(center, ncol = 1))) + 
  kronecker(matrix(1, 1, M), matrix(N_shift, ncol = 1))
