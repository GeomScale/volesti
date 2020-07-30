#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope


print("I imported everything right!")

if __name__ == "__main__":
    dim = 8
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


# Run tests for the compute_volume() function
    volume_SoB = p.compute_volume(vol_method="sequence_of_balls", walk_method="uniform_ball", walk_len=5, epsilon=0.05, seed=volestipy.get_time_seed())
    print("test1.py: Volume (sequence of balls): {}".format(volume_SoB))

    volume_GA  = p.compute_volume(vol_method = "cooling_gaussian", walk_method="gaussian_CDHR", walk_len=5, epsilon=0.05, seed=42)
    print("test1.py: Volume (gaussian annealing): {}".format(volume_GA))

    test_volume = p.compute_volume()
    print("test1.py: Volume (test_volume): {}".format(test_volume))

# Run tests for the generate_samples() function
    samples = p.generate_samples(walk_len = 1, number_of_points = 1000, number_of_points_to_burn = 50, boundary = True, cdhr = True, rdhr = False, 
     gaussian = False, set_L = False, billiard = False, ball_walk = False, a = 0, L = 0)
    print("test1.py: Samples:")
    print(samples)


    rounding = p.rounding(walk_len = 2, billiard = True)
    print("rounding function ran ok\n")
    print(rounding[2])



