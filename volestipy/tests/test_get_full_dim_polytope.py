#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope
from volestipy import low_dim_HPolytope

if __name__ == "__main__":
    
    m = 3
    n = 5

    A = np.zeros((2*n, n), dtype=np.float)
    A[0:n] = np.eye(n)
    A[n:] -=  np.eye(n,n, dtype=np.float)
    print("\n This is the A matrix: \n")
    print(A)

    b = np.ones(2*n, dtype=np.float)
    print("\n This is the vector b: \n")
    print(b)

    Aeq = np.random.normal(0, 1, size=(m, n))
    print("\n This is the Aeq matrix: \n")
    print(Aeq)

    beq = np.zeros(m)
    print("\n This is the vector beq: \n")
    print(beq)

    # run the low_dim_HPolytope() class and build an object 
    low_hp = low_dim_HPolytope(A,b,Aeq,beq)
    # run the full_dimensiolal_polytope() function
    get_full_hp = low_hp.full_dimensiolal_polytope()

    # print the output of the full_dimensiolal_polytope() function
    print("\n this is the get_full_hp.A, meaning the A matrix of the full dimensional polytope:")
    print(get_full_hp.A)
    print("\n this is the get_full_hp.b, meaning the b vector of the full dimensional polytope:")
    print(get_full_hp.b)
