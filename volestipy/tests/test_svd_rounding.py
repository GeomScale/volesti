#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope

if __name__ == "__main__":
    dim = 5
    A = np.zeros((2*dim, dim), dtype=np.float)
    A[0:dim] = np.eye(dim)
    A[dim:] -=  np.eye(dim,dim, dtype=np.float)
    b = np.ones(2*dim, dtype=np.float)


    hp = HPolytope(A,b)

#    print("The polytope:")
#    print("Dimensions: {}".format(hp.dimensions))
#    print("A matrix:")
#    print(hp.A)
#    print(type(hp.A))
#    print("b vector:")
#    print(hp.b)
#    print(type(hp.b))

    # Run test for the rounding() function and its different methods
    rounding_returns = ["new_A","new_b","T_matrix","shift","round_val"]
#    print("*** This is a test for th rounding() function of the volestipy library *** \n")

    # Case 1
    rounding = hp.rounding_svd()
#    print("\n >> This is the output for the svd rounding method <<")
#    print(rounding_output_svd[0])


    print(rounding[0])
