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
    print("*** This is a test for th compute_volume() function of the volestipy library *** \n")
    print("The polytope:")
    print("Dimensions: {}".format(p.dimensions))
    print("A matrix:")
    print(p.A)
    print("b vector:")
    print(p.b)

# Run tests for the compute_volume() function and its different methods
    print("\n>> Here are the outcomes of the different methods for computing volume <<\n")
    volume_SoB = p.compute_volume(vol_method = "sequence_of_balls", walk_method = "uniform_ball", walk_len = 5, epsilon = 0.05, seed = volestipy.get_time_seed())
    print("Volume (sequence of balls): {}\n".format(volume_SoB))

    volume_GA  = p.compute_volume(vol_method = "cooling_gaussian", walk_method = "gaussian_CDHR", walk_len = 5, epsilon = 0.05, seed = 42)
    print("Volume (gaussian annealing): {} \n".format(volume_GA))

    test_volume = p.compute_volume(vol_method = "cooling_balls", walk_method = "uniform_ball", walk_len = 5, epsilon = 0.05, seed = 42)
    print("Volume (test_volume): {}\n".format(test_volume)) 


