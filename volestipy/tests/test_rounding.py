#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope

if __name__ == "__main__":
    dim = 9
    A = np.zeros((2*dim, dim), dtype=np.float)
    A[0:dim] = np.eye(dim)
    A[dim:] -=  np.eye(dim,dim, dtype=np.float)
    b = np.ones(2*dim, dtype=np.float)


    hp = HPolytope(A,b)

    print("The polytope:")
    print("Dimensions: {}".format(hp.dimensions))
    print("A matrix:")
    print(hp.A)
    print(type(hp.A))
    print("b vector:")
    print(hp.b)
    print(type(hp.b))

    # Run test for the rounding() function and its different methods
    rounding_returns = ["new_A","new_b","T_matrix","shift","round_val"]

    # Case 1
    rounding_output_max_ellipsoid = hp.rounding(rounding_method = "max_ellipsoid")
    print("\n ***This is the output for the max_ellipsoid rounding method.***\n")
    for i in range(len(rounding_output_max_ellipsoid)):
        print("\n" + rounding_returns[i] + ": ")
        print(rounding_output_max_ellipsoid[i])

    # Case 2
    rounding_output_svd = hp.rounding(rounding_method = "svd")
    print("\n ***This is the output for the svd rounding method.***\n")
    for i in range(len(rounding_output_svd)):
        print("\n" + rounding_returns[i] + ": ")
        print(rounding_output_svd[i])


    # Case 3
    rounding_output_min_ellipsoid = hp.rounding(rounding_method = "min_ellipsoid")
    print("\n ***This is the output for the min_ellipsoid rounding method.***\n")
    for i in range(len(rounding_output_min_ellipsoid)):
        print("\n" + rounding_returns[i] + ": ")
        print(rounding_output_min_ellipsoid[i])

