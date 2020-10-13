#!/usr/bin/python3

import numpy as np
import gurobipy as gp
from volestipy import *
import matplotlib.pyplot as plt
import sys, datetime

start = datetime.datetime.now()

# Set a variable with the input / metabolic network file
input_file = '../../bigg_files/iYO844.json'

# Read json
read_ecoli_core = read_json_file(input_file)

# Pre-process it
A = read_ecoli_core[0]
b = read_ecoli_core[1]
Aeq = read_ecoli_core[2]
beq = read_ecoli_core[3]

# Pre-process it
proc = pre_process(A, b, Aeq, beq)
A_proc = proc[0]
b_proc = proc[1]
Aeq_proc = proc[2]
beq_proc = proc[3]
min_fluxes = proc[4]
max_fluxes = proc[5]

# Get an object for the low_dim_HPolytope class for the pre-processed polytope
low_hp = low_dim_HPolytope(A_proc, b_proc, Aeq_proc, beq_proc)

## And then get the full dimensional polytope
get_fd_hp = low_hp.full_dimensiolal_polytope()
A_fd = get_fd_hp[0].A
b_fd = get_fd_hp[0].b
N = get_fd_hp[1]
N_shift = get_fd_hp[2]

# Get the max ball for the full dimensional polytope
max_ball_center_point, max_ball_radius = get_max_ball(A_fd, b_fd)
print("max ball was calculated for the full dimensional polytope.")

### Now we can use the full dimensional polytope; but before sampling on it, we need to round it

# First, initialize an HPolytope using the full dimensional polytope features we got
hp = HPolytope(A_fd, b_fd)

## Then use one of the volestipy functions for rounding
rounding_returns = ["new_A","new_b","T_matrix","shift","round_val"]

# Rounding by making use of max ball and the max_ellipsoid method
rounding_output_svd = hp.rounding_svd()
rounded_A = rounding_output_svd[0]
rounded_b = rounding_output_svd[1]
rounded_T = rounding_output_svd[2]
rounded_shift = rounding_output_svd[3]
print("Rounding has been completed")

## Finally, generate random samples from the rounded full dimensional polytope

# Get the max ball for the rounded polytope
rounded_center_point, rounded_radius = get_max_ball(rounded_A, rounded_b)

# Build the rounded polytope
rounded_polytope = HPolytope(rounded_A, rounded_b)

# and calculate L parameter for sampling
d = rounded_polytope.dimensions
L_value = 4 * d * rounded_radius

# Then, sample on it
samples = rounded_polytope.generate_samples(walk_len = 5, number_of_points = 10000, number_of_points_to_burn = 50, \
                                            radius = rounded_radius, inner_point = rounded_center_point, L = L_value)

# Now map the points retrieved to the initial polytope
mapped_samples = map_samples_on_initial_polytope(samples, rounded_T, rounded_shift, N, N_shift)

# Plot the distribution of the flux unde study, e.g ATPS4r
fluxes = mapped_samples[21,:]
plt.hist(fluxes, bins='auto', histtype='bar', rwidth=0.7)
plt.savefig('flux_22.png')

end = datetime.datetime.now()
total_time = end - start
print("Script totally ran for :", total_time)
