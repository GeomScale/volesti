#!/usr/bin/python3.6

import numpy as np
import gurobipy as gp
from volestipy import *
import sys

input_file = 'e_coli_core.json'


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


print(A_fd)
print(b_fd)


### Now we can use the full dimensional polytope; but before sampling on it, we need to round it

# First, initialize an HPolytope using the full dimensional polytope features we got
hp = HPolytope(A_fd, b_fd)

## Then use one of the volestipy functions for rounding
#rounding_output_svd = hp.rounding(rounding_method = "svd")
#
#print("This is the shape of the full dimensional A: ")
#print(A_fd[0].A.shape)



