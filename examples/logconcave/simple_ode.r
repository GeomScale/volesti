# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou

# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for ODE solvers 
library(volesti)

F <- function (x) (-x)
order <- 2
step_size <- 0.01
n <- 1000
initial_conditions <- list("x_1" = c(0), "x_2" = c(1))
initial_time <- 0

# Do not impose constraints
domains <- list()

# Call the ode solver
states <- volesti::ode_solve(dimension=1, n=n, F=F, initial_time=initial_time, step_size=step_size, order=order, method="leapfrog", initial_conditions=initial_conditions, domains = list())

x <- states[["x_1"]]
v <- states[["x_2"]]

plot(x, v)

