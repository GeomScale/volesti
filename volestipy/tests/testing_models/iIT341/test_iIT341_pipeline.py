#!/usr/bin/python3

import numpy as np
import gurobipy as gp
from volestipy import *
import matplotlib.pyplot as plt
import sys, datetime
from scipy import linalg 

np.set_printoptions(threshold=sys.maxsize)

start = datetime.datetime.now()

# Set a variable with the input / metabolic network file
input_file = '../../bigg_files/iSDY_1059.json'
print("input_file is: ") ; print(input_file)


# Read json
read_ecoli_core = read_json_file(input_file)

# Pre-process it
A = read_ecoli_core[0]
b = read_ecoli_core[1]
Aeq = read_ecoli_core[2]
beq = read_ecoli_core[3]

# Pre-process it
proc = pre_process(A, b, Aeq, beq)
A_proc = proc[0] ; b_proc = proc[1] ; Aeq_proc = proc[2] ; beq_proc = proc[3] ; min_fluxes = proc[4] ; max_fluxes = proc[5]

N = linalg.null_space(Aeq_proc)
N_shift = np.linalg.lstsq(Aeq_proc, beq_proc, rcond=None)[0]

product = np.dot(A_proc, N_shift)
b_fd = np.subtract(b_proc, product)

A_fd = np.dot(A_proc,N)
lines = []
b_elements = []

for i in range(A_fd.shape[0]):
      entry = np.linalg.norm(A_fd[i,])

      if entry < 1e-06:
         continue
      else:
         lines.append(A_fd[i,:])
         b_elements.append(b_fd[i])

A_fd_true = np.array(lines)
b_fd_true = np.array(b_elements)

np.save('recon_A_fd_true.npy', A_fd_true)
np.save('recon_b_fd_true.npy', b_fd_true)
np.save('recon_n_shift.npy', N_shift)

# Build an object of the full dimensional polytope
hp_true = HPolytope(A_fd_true, b_fd_true)

####        FIRST APPROACH

try:
   b_fd = b_proc
   # Build an object of the full dimensional polytope
   hp = HPolytope(A_fd, b_fd)
   
   # Get the max ball for the full dimensional polytope
   max_ball_center_point, max_ball_radius = get_max_ball(A_fd, b_fd)
   print("max ball center pointer for NON-scaled polytope before rounding is: ") ; print(max_ball_center_point)
   print("max ball radius for NON-scaled polytope  before rounding is: ") ; print(max_ball_radius)

except:
   print("Cannot get max ball with b_fd = b_proc") 

print("\n\n\n-----------------------------------------------------------\n\n\n\n")

####        SECOND APPROACH

try:
  
   # Get the max ball for the full dimensional polytope
   max_ball_center_point, max_ball_radius = get_max_ball(A_fd_true, b_fd_true)
   print("max ball center pointer for NON-scaled TRUE polytope before rounding is: ") ; print(max_ball_center_point)
   print("max ball radius for NON-scaled TRUE polytope before rounding is: ") ; print(max_ball_radius)
   print(A_fd.shape[1]) 
   
except Exception:
   print("Cannot get max ball with  b_proc - product where product = np.dot(A_proc, N_shift) ") 

print("\n\n\n-----------------------------------------------------------\n\n\n\n")

####        THIRD APPROACH
   
# Check whether the max_ball_center_point is a zero-array and if yes it is try with scaling it

# if max_ball_radius <= 0:

try:
    
   # And scale it
   cscale, rscale = gmscale(A_fd_true, 5, 0.90)
   scaled_A, scaled_b, diag_matrix = scaled_polytope(hp_true, cscale, rscale)
   
   # Now build a new object for the scaled full dimensional polytope
   hp_scaled = HPolytope(scaled_A, scaled_b)
   
   # Get the max ball for the full dimensional polytope
   print("Computing max ball...")
   scaled_max_ball_center_point, scaled_max_ball_radius = get_max_ball(scaled_A, scaled_b)
   
   print("max ball center pointer for scaled TRUE polytope before rounding is: ") ; print(scaled_max_ball_center_point)
   print("max ball radius for scaled TRUE polytope before rounding is: ") ; print(scaled_max_ball_radius)

except:
   print("Sorry. I cannot deal with this metabolic network.")

sys.exit(0)










## Then use one of the volestipy functions for rounding
rounding_returns = ["new_A","new_b","T_matrix","shift","round_val"]

# Rounding by making use of max ball and the max_ellipsoid method
print("Rounding is about to start")
rounding_output_svd = hp_scaled.rounding_svd(scale = diag_matrix)
rounding_output_svd = hp.rounding_svd()
#rounding_output_svd = hp.rounding(rounding_method = 'max_ellipsoid')

rounded_A = rounding_output_svd[0]
rounded_b = rounding_output_svd[1]
rounded_T = rounding_output_svd[2]
rounded_shift = rounding_output_svd[3]
print("rounding has been completed")

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
