#!/usr/bin/python3

import numpy as np
import gurobipy as gp
from volestipy import *

# Set a variable with the input / metabolic network file
input_file = 'bigg_files/e_coli_core.json'

# Read json
read_ecoli_core = read_json_file(input_file)

# Pre-process it
A = read_ecoli_core[0]
b = read_ecoli_core[1]
Aeq = read_ecoli_core[2]
beq = read_ecoli_core[3]

## Pre-process it
#proc = pre_process(A, b, Aeq, beq)
#A_proc = proc[0]
#b_proc = proc[1]
#Aeq_proc = proc[2]
#beq_proc = proc[3]


#A_proc = np.load('A_preprocessed.npy')
#b_proc = np.load('b_preprocessed.npy')
#Aeq_proc = np.load('Aeq_preprocessed.npy')
#beq_proc = np.load('beq_preprocessed.npy')




# Get an object for the low_dim_HPolytope class for the pre-processed polytope
low_hp = low_dim_HPolytope(A, b, Aeq, beq)
print("object ok")


# And then get the full dimensional polytope
get_fd_hp = low_hp.full_dimensiolal_polytope()
A_fd = get_fd_hp[0].A
b_fd = get_fd_hp[0].b
N = get_fd_hp[1]
N_shift = get_fd_hp[2]

print("\n\n *** This is the full dimensional polytope ***")
print(A_fd)
print(b_fd)

# Get the max ball for the full dimensional polytope
print("\n\n\n We are about to calculate max ball")
max_ball_center_point, max_ball_radius = get_max_ball(A_fd, b_fd)
print("Max ball was found:")
print(max_ball_center_point)
print(max_ball_radius)

### Now we can use the full dimensional polytope; but before sampling on it, we need to round it

# First, initialize an HPolytope using the full dimensional polytope features we got
hp = HPolytope(A_fd, b_fd)
print("An HP was built out of the full dimensional polytope features")

## Then use one of the volestipy functions for rounding
rounding_returns = ["new_A","new_b","T_matrix","shift","round_val"]


# Rounding by making use of max ball and the max_ellipsoid method
print("Rounding with max ball and max_ellipsoid \n")
rounding_output_max_ellipsoid = hp.rounding(rounding_method = "max_ellipsoid", inner_point = max_ball_center_point, radius = max_ball_radius)
rounded_A = rounding_output_max_ellipsoid[0]
rounded_b = rounding_output_max_ellipsoid[1]
rounded_T = rounding_output_max_ellipsoid[2]
rounded_shift = rounding_output_max_ellipsoid[3]

for i in range(len(rounding_output_max_ellipsoid)):
   print("\n" + rounding_returns[i] + ": ")
   print(rounding_output_max_ellipsoid[i])

## Finally, generate random samples from the rounded full dimensional polytope

# Get the max ball for the rounded polytope
rounded_center_point, rounded_radius = get_max_ball(rounded_A, rounded_b)

# Build the rounded polytope
rounded_polytope = HPolytope(rounded_A, rounded_b)


# Then, sample on it
samples = rounded_polytope.generate_samples(walk_len = 5, number_of_points = 10000, number_of_points_to_burn = 50, radius = rounded_radius, inner_point = rounded_center_point)
print("\n >> This is the output for random sampling algorithm using max ball on a rounded polytope <<\n")
print("Samples:")
print(samples)
print(type(samples))


# Now map the points retrieved to the initial polytope
print("Points sampled are under mapping")
final_output = map_samples_on_initial_polytope(samples, rounded_T, rounded_shift, N, N_shift)
print("These are the final points")
print(final_output)

