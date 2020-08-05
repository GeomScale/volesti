#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope

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





