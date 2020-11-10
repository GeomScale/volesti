#!/usr/bin/python3

import numpy as np
import gurobipy as gp
from volestipy import *
import sys, os
from os import listdir
from os.path import isfile, join
from scipy import linalg 


# Make a list with all the BIGG model files
mypath = '/home/tolis/data/metabolic_json/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
#onlyfiles = '/home/tolis/data/metabolic_json/e_coli_core.json'


def run_pipeline(input_file):

   name_of_met_net = input_file[:-5]
   output_dir = os.getcwd()
   output_dir = output_dir + '/volestipy_output/'
   output_dir_net = output_dir + name_of_met_net
   
   # Check if directories already made   
   if os.path.isdir(output_dir) == False:
      os.mkdir(output_dir, exist_ok=True)

   if os.path.isdir(output_dir_net) == False:
      os.mkdir(output_dir_net)
   
   # Move to the model's output directory
   os.chdir(output_dir_net)
   
   # Read json
   input_file = mypath + input_file
   met_net = read_json_file(input_file)
   
   # Keep the initial data from the model
   A = met_net[0] ; b = met_net[1] ; Aeq = met_net[2] ; Aeq = np.ascontiguousarray(Aeq) ; beq = met_net[3]
   
   # Pre-process it
   proc = pre_process(A, b, Aeq, beq)
   A_proc = proc[0] ; b_proc = proc[1] ; Aeq_proc = proc[2] ; beq_proc = proc[3]
   

   # [ATTENTION!] New part on the pipeline - replacing the get full dimensional function
   
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

   np.save(os.path.join(output_dir_net,'A_fd_true.npy'), A_fd_true)
   np.save(os.path.join(output_dir_net,'b_fd_true.npy'), b_fd_true)
   np.save(os.path.join(output_dir_net,'n_shift.npy'), N_shift)
   np.save(os.path.join(output_dir_net,'n.npy'), N)   
   
   # Build an object of the full dimensional polytope
   hp_true = HPolytope(A_fd_true, b_fd_true)
      
   
   ## -------- 3 approaches to get max ball ------
   
   approach_1 = False ; approach_2 = True ; approach_3 = True # approach_1 is False in terms of not having the approach 1 running in this test
   
   ####        FIRST APPROACH - no scale , no true polytope
   
   # try:
   #    low_hp = low_dim_HPolytope(A, b, Aeq, beq)
   #    get_fd_hp = low_hp.full_dimensiolal_polytope()
   #    A_fd = get_fd_hp[0].A
   #    b_fd = get_fd_hp[0].b
   #    N = get_fd_hp[1]
   #    N_shift = get_fd_hp[2]
   #    
   #    # Make b full dimensional equal to the processed one
   #    b_fd = b_proc
   #    
   #    # Get the max ball for the full dimensional polytope
   #    max_ball_center_point, max_ball_radius = get_max_ball(A_fd, b_fd)
   # 
   #    # print("max ball center pointer for NON-scaled polytope before rounding is: ") ; print(max_ball_center_point)      
   #    print("Approach 1: no scale, no true polytope, b_fd = b_proc")
   #    print("max ball radius for NON-scaled polytope, NO true polytope, with b_fd = b_proc: ") ; print(max_ball_radius)
   # 
   # except:
   #    approach_1 = False
   #    print("Cannot get max ball with b_fd = b_proc") 
   # 
   # print("\n\n\n-----------------------------------------------------------\n\n\n\n")
   
   ####        SECOND APPROACH - no scale, true polytope
   
   try:
   
      # Get the max ball for the full dimensional polytope
      max_ball_center_point, max_ball_radius = get_max_ball(A_fd_true, b_fd_true)
   
      # Test if we have an error that does not make our program to fail
      b_check = b_fd_true - np.dot(A_fd_true, max_ball_center_point)
      product = ( 1 / ( (1 / max_ball_radius ) ** ( 1 /A_fd_true.shape[1] )))
   
      A_check = A_fd_true * ( 1 / ( (1 / max_ball_radius ) ** ( 1 / A_fd_true.shape[1] )))
      check_center, check_radius = get_max_ball(A_check, b_check)
   
      A_check_half = A_fd_true * 0.5
      check_center_half, check_radius_half = get_max_ball(A_check, b_check)
      
      # print("\n\nmax ball center pointer for NON-scaled TRUE polytope before rounding is: ") ; print(max_ball_center_point)
      print("Approach 2: no scale, true polytope")
      print("product is: ") ; print(product)     
      print("check_radius") ; print(check_radius)
      print("check_half_radius: ") ; print(check_radius_half)
      
      print("Actual radius returned: ") ; print(max_ball_radius)
   
   
   except Exception:
      approach_2 = False
      print("Cannot get max ball with  b_proc - product where product = np.dot(A_proc, N_shift) ") 
   
   print("\n\n\n-----------------------------------------------------------\n\n\n\n")
   

   ####        THIRD APPROACH - scale and true 

   try:
   
      # And scale it
      cscale, rscale = gmscale(A_fd_true, 5, 0.90)
      scaled_A, scaled_b, diag_matrix = scaled_polytope(hp_true, cscale, rscale)
   
      # Now build a new object for the scaled full dimensional polytope
      hp_scaled = HPolytope(scaled_A, scaled_b)
   
      # Get the max ball for the full dimensional polytope
      scaled_max_ball_center_point, scaled_max_ball_radius = get_max_ball(scaled_A, scaled_b)

      # Test if we have an error that does not make our program to fail
      b_check = scaled_b - np.dot(scaled_A, scaled_max_ball_center_point)
      product_2 = ( 1 / ( (1 / scaled_max_ball_radius ) ** ( 1 / scaled_A.shape[1] )))
   
      A_check = scaled_A * product_2
      check_center_scaled, check_radius_scaled = get_max_ball(A_check, b_check)
   
      A_check_half = scaled_A * 0.5
      check_center_scaled_half, check_radius_scaled_half = get_max_ball(A_check, b_check)
      
      
      print(scaled_max_ball_center_point)
      print("Approach 3: scaled and true")
      print("product is: ") ; print(product)
      print("check_radius") ; print(check_radius_scaled)
      print("check_radius_scaled_half: ") ; print(check_radius_scaled_half)

      print("Actual radius returned: ") ; print(max_ball_radius)
   
   except:
      approach_3 = False
      print("Sorry. I cannot deal with this metabolic network.")
   

   ####
   ####  Investigating for max_ball has been completed !!! move to rounding step. 
   ####

   # Rounding using the greatest approach

   if approach_3 == True:
      
      rounding_svd_output = hp_scaled.rounding_svd(scale = diag_matrix)
      rounded_A = rounding_svd_output[0] ; rounded_b = rounding_svd_output[1] ; rounded_T = rounding_svd_output[2] ; rounded_shift = rounding_svd_output[3]
      
      rounded_shift = rounded_shift + scaled_max_ball_center_point
      rounded_T = rounded_T * product_2
     
      
      
      print("model " + name_of_met_net + "SVD rounding completed with approach 3. \n")

   elif approach_2 == True and approach_1 == False: 

      rounding_svd_output = hp_true.rounding_svd()
      rounded_A = rounding_svd_output[0]
      rounded_b = rounding_svd_output[1]
      rounded_T = rounding_svd_output[2]
      rounded_shift = rounding_svd_output[3]
      print("model " + name_of_met_net +  "SVD rounding completed with approach 2. \n")
      
   elif approach_1 == True:
      
      hp = HPolytope(A_fd, b_fd)
      
      roundign_svd_output = hp.rounding_svd()
      rounded_A = rounding_svd_output[0]
      rounded_b = rounding_svd_output[1]
      rounded_T = rounding_svd_output[2]
      rounded_shift = rounding_svd_output[3]
      print("model " + name_of_met_net + "SVD rounding completed with approach 1. \n")
   
   
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
   
   sys.exit(0)

# --------------------------------------------------------------------------------

# Run the pipeline at last! 
for model in onlyfiles:
   run_pipeline(model)











   ##  ------ this has been replaced  -------
   
   # # Get an object for the low_dim_HPolytope class for the pre-processed polytope
   # low_hp = low_dim_HPolytope(A, b, Aeq, beq)
   # print("object ok")
   
   ## And then get the full dimensional polytope
   #get_fd_hp = low_hp.full_dimensiolal_polytope()
   #A_fd = get_fd_hp[0].A
   #b_fd = get_fd_hp[0].b
   #N = get_fd_hp[1]
   #N_shift = get_fd_hp[2]
   
   #print("\n\n *** This is the full dimensional polytope ***")
   #print(A_fd)
   #print(b_fd)
   
   # A_fd = np.load('A_full_dim.npy')
   # b_fd = np.load('b_full_dim.npy')
   # N = np.load('N_full_dim.npy')
   # N_shift = np.load('shift_full_dim.npy')
   
   # # Get the max ball for the full dimensional polytope
   # print("\n\n\n We are about to calculate max ball")
   # max_ball_center_point, max_ball_radius = get_max_ball(A_fd, b_fd)
   # print("Max ball was found:")
   # print(max_ball_center_point)
   # print(max_ball_radius)
   
   ##  ------ this has been replaced  -------


