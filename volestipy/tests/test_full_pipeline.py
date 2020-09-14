#!/usr/bin/python3.6

import numpy as np
import gurobipy as gp
from volestipy import *
import sys

input_file = 'bigg_files/e_coli_core.json'


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

# Get an object for the low_dim_HPolytope class for the pre-processed polytope
low_hp = low_dim_HPolytope(A_proc, b_proc, Aeq_proc, beq_proc)

# And then get the full dimensional polytope
get_fd_hp = low_hp.full_dimensiolal_polytope()
A_fd = get_fd_hp[0].A
b_fd = get_fd_hp[0].b

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
print("\n\n Rounding is about to start")
rounding_output_min_ellipsoid = hp.rounding(rounding_method = "min_ellipsoid", inner_point = max_ball_center_point, radius = max_ball_radius)

rounding_returns = ["new_A","new_b","T_matrix","shift","round_val"]
for i in range(len(rounding_output_min_ellipsoid)):
   print("\n" + rounding_returns[i] + ": ")
   print(rounding_output_min_ellipsoid[i])

# Check for the rest rounding methods
rounding_output_max_ellipsoid = hp.rounding(rounding_method = "max_ellipsoid")
#for i in range(len(rounding_output_max_ellipsoid)):
#   print("\n" + rounding_returns[i] + ": ")
#   print(rounding_output_max_ellipsoid[i])


