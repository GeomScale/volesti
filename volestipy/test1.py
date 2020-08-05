#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope

if __name__ == "__main__":
    dim = 3
    A = np.zeros((2*dim, dim), dtype=np.float)
    A[0:dim] = np.eye(dim)
    A[dim:] -=  np.eye(dim,dim, dtype=np.float)
    b = np.ones(2*dim, dtype=np.float)


    p = HPolytope(A,b)
    print("The polytope:")
    print("Dimensions: {}".format(p.dimensions))
    print("A matrix:")
    print(p.A)
    print(type(p.A))
    print("b vector:")
    print(p.b)
    print(type(p.b))


# Run tests for the compute_volume() function and its different methods
    volume_SoB = p.compute_volume(vol_method="sequence_of_balls", walk_method="uniform_ball", walk_len=5, epsilon=0.05, seed=volestipy.get_time_seed())
    print("test1.py: Volume (sequence of balls): {}".format(volume_SoB))

    volume_GA  = p.compute_volume(vol_method = "cooling_gaussian", walk_method="gaussian_CDHR", walk_len=5, epsilon=0.05, seed=42)
    print("test1.py: Volume (gaussian annealing): {}".format(volume_GA))

    test_volume = p.compute_volume()
    print("test1.py: Volume (test_volume): {}".format(test_volume))
    print(" \n\n ***The computing of the volume of the polytope has been concluded. *** \n\n")
    print("---------------------------------------------------------")

# Run tests for the generate_samples() function
    samples = p.generate_samples(walk_len = 1, number_of_points = 1000, number_of_points_to_burn = 50, boundary = True, cdhr = True, rdhr = False,
     gaussian = False, set_L = False, billiard = False, ball_walk = False, a = 0, L = 0)
    print("test1.py: Samples:")
    print(samples)
    print("\n\n ***The sampling step is over. ***\n\n")
    print("----------------------------------------------------------")

# Run test for the rounding() function and its different methods
    print("\n\n ***The rounding step is about to start*** \n\n")
    p.rounding(rounding_method = "max_ellipsoid")
    print("\n ***This is the output for the max_ellipsoid rounding method.***\n")



#    for i in rounding_output_max_ellipsoid:
#        print(i)
#        print("\n\n*************\n\n")
#
#    rounding_output_svd = p.rounding(rounding_method = "svd")
#    print("\n\nthis is the output for the svd rounding")
#    for i in rounding_output_svd:
#        print(i)
#        print("\n\n*************\n\n")
#
#
#    rounding_output_min_ellipsoid = p.rounding(rounding_method = "min_ellipsoid")
#    print("\n this is the output for the max_ellipsoid rounding method\n")
#    for i in rounding_output_min_ellipsoid:
#        print(i)
#        print("\n\n*************\n\n")
