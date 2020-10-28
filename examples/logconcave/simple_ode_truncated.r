# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou

# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for truncated ODE solvers 
library(volesti)

F <- function (x) (x)
order <- 1
step_size <- 0.01
n <- 1000
initial_conditions <- list("x_1" = c(0.1))
initial_time <- 0

A <- matrix(c(1, -1), ncol=1, nrow=2, byrow=TRUE)
b <- c(1, 0)

# Create domain of truncation
P_1 <- volesti::Hpolytope$new(A, b)
domains <- list("P_1" = P_1)

# Call the ode solver
states <- volesti::ode_solve(dimension=1, n=n, F=F, initial_time=initial_time, step_size=step_size, order=order, method="euler", initial_conditions=initial_conditions, domains = domains)

x <- states[["x_1"]]
t <- step_size * seq(0, n - 1)

plot(t, x)
