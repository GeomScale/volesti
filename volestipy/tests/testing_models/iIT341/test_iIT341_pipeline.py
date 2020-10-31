#!/usr/bin/python3

import numpy as np
import gurobipy as gp
from volestipy import *
import matplotlib.pyplot as plt
import sys, datetime
from scipy import linalg 


# np.set_printoptions(threshold=sys.maxsize)

start = datetime.datetime.now()

# Set a variable with the input / metabolic network file
input_file = '../../bigg_files/iIT341.json'

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

# print("Aeq_proc shape: " ) ; print(Aeq_proc.shape)
# print("beq_proc shape: ") ; print(beq_proc.shape)
# print("b_proc type is: ") ; print(type(b_proc))

## Get an object for the low_dim_HPolytope class for the pre-processed polytope
# low_hp = low_dim_HPolytope(A_proc, b_proc, Aeq_proc, beq_proc)

## And then get the full dimensional polytope
# get_fd_hp = low_hp.full_dimensiolal_polytope()


# N = get_fd_hp[1]
N = linalg.null_space(Aeq_proc)
print("N shape is :") ; print(N.shape)


# N_shift = get_fd_hp[2]
# N_shift = np.linalg.solve(Aeq_proc, beq_proc)
# N_shift = np.zeros(Aeq_proc.shape[1])
N_shift = np.linalg.lstsq(Aeq_proc, beq_proc, rcond=None)[0]
# print("N_shift is :") ; print(type(N_shift))



# A_fd = get_fd_hp[0].A
A_fd = np.dot(A_proc,N)
# print("A_fd shape is :") ; print(A_fd.shape)

# b_fd = get_fd_hp[0].b


# print("**************************************************")
# print("A_proc.shape: ") ; print(A_proc.shape)
# print("**************************************************")
# print("N_shift type: ") ; print(type(N_shift))
# print("N_shift.shape: ") ; print(N_shift.shape)
# # print("N_shift: ") ; print(N_shift)
# print("**************************************************")
# print("b_proc.shape: ") ; print(b_proc.shape)
# print("b_proc is: ") ; print(b_proc)
# print("**************************************************")
# print("beq_proc type: ") ; print(type(beq_proc))
# print("beq_proc shape: ") ; print(beq_proc.shape)
# # print("beq_proc: ") ; print(beq_proc)
# print("**************************************************")

try:
   b_fd = b_proc
   # Build an object of the full dimensional polytope
   hp = HPolytope(A_fd, b_fd)
   
   # Get the max ball for the full dimensional polytope
   max_ball_center_point, max_ball_radius = get_max_ball(A_fd, b_fd)
   # print("max ball center pointer for NON-scaled polytope before rounding is: ") ; print(max_ball_center_point)
   # print("max ball radius for NON-scaled polytope  before rounding is: ") ; print(max_ball_radius)

except:
   print("Cannot get max ball with b_fd = b_proc") 

try:

   N_shift = np.zeros(N_shift.shape[0])
   print("N_shift: ") ; print(N_shift)

   product = np.dot(A_proc, N_shift)
   print("product A_proc * N_shift : ") ; print(product)

   print("b_proc : ") ; print(b_proc)

   b_fd = np.subtract(b_proc, product)
   print("b_fd shape is with product is:") ; print(b_fd.shape)
   print("b_fd is with product is:") ; print(b_fd)   
   

   print("A_fd") ; print(A_fd)
   A_fd = 0.7 * A_fd
   print("multiplied A_fd") ; print(A_fd)

   # Build an object of the full dimensional polytope
   hp = HPolytope(A_fd, b_fd)
   
   # Get the max ball for the full dimensional polytope
   max_ball_center_point, max_ball_radius = get_max_ball(A_fd, b_fd)
   # print("max ball center pointer for NON-scaled polytope before rounding is: ") ; print(max_ball_center_point)
   # print("max ball radius for NON-scaled polytope  before rounding is: ") ; print(max_ball_radius)

except Exception:
   print("Cannot get max ball with  b_proc - product where product = np.dot(A_proc, N_shift) ") 



sys.exit(0)



   
# Check whether the max_ball_center_point is a zero-array and if yes it is try with scaling it
if max_ball_radius <= 0:

   try:
       
      # And scale it
      cscale, rscale = gmscale(A_fd, 5, 0.99)
      scaled_A, scaled_b, diag_matrix = scaled_polytope(hp, cscale, rscale)
      
      print("scaled_A: ") ; print(scaled_A)
      print("scaled_b: ") ; print(scaled_b)
      
      # Now build a new object for the scaled full dimensional polytope
      hp_scaled = HPolytope(scaled_A, scaled_b)
      
      # Get the max ball for the full dimensional polytope
      print("Computing max ball...")
      scaled_max_ball_center_point, scaled_max_ball_radius = get_max_ball(scaled_A, scaled_b)
      
      
      print("max ball center pointer for scaled polytope before rounding is: ") ; print(scaled_max_ball_center_point)
      print("max ball radius for scaled polytope  before rounding is: ") ; print(scaled_max_ball_radius)

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
