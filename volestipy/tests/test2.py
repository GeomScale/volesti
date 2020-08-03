#!/usr/bin/python3.6

import numpy as np
from volestipy import HPolytope


if __name__ == "__main__":
    # cube4d
    # b = [1, 0, 1, 0, 1, 0, 1, 0]
    # A =[[1, 0, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, -1, 0], [0, 0, 0, 1], [0, 0, 0, -1]]

    # #cubetriangles4d
    b = [1, 0, 1, 0, 1, 0, 1, 0, 0]
    A =[[1, 0, 0, 0], [-1, 0, 0, 0], [0, 1, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, -1, 0], [0, 0, 0, 1], [0, 0, 0, -1], [1, -1, 0, 0]]

    #cube2d
    # b = [1, 1, 1, 1]
    # A =[[1, 0], [-1, 0], [0, 1], [0, -1]]


    # cube#trinagles2d
    # b = [1, 1, 1, 1, -1]
    # A =[[1, 0], [-1, 0], [0, 1], [0, -1], [1,-1]]

    A = np.array(A, dtype=np.float)
    b = np.array(b, dtype=np.float)


    p = HPolytope(A,b)
    print("The polytope:")
    print("Dimensions: {}".format(p.dimensions))
    print("A matrix:")
    print(p.A)
    print("b vector:")
    print(p.b)

    volume_SoB = p.compute_volume(walk_len=5, epsilon=0.05, vol_method="sequence_of_balls", seed=42)
    volume_GA = p.compute_volume(walk_len=5, epsilon=0.05, vol_method = "cooling_gaussian", seed=42)
    samples = p.generate_samples(walk_len = 5, number_of_points = 80000, number_of_points_to_burn = 50, boundary = True, cdhr = True, rdhr = False, 
              gaussian = False, set_L = False, billiard = False, ball_walk = False, a = 0, L = 0)



    print("Volume (sequence of balls): {}".format(volume_SoB))
    print("Volume (gaussian annealing): {}".format(volume_GA))
    print("Samples:")
    print(samples)
