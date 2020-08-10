#!/usr/bin/python3.6

import numpy as np
import volestipy
from volestipy import HPolytope

if __name__ == "__main__":
    
    m = 3
    n = 5

    A = np.zeros((2*n, n), dtype=np.float)
    A[0:n] = np.eye(n)
    A[n:] -=  np.eye(n,n, dtype=np.float)
    print("This is the A matrix: \n")
    print(A)

    b = np.ones(2*n, dtype=np.float)
    print("This is the vector b: \n")
    print(b)

    Aeq = np.random.normal(1, size=(m, n))
    print("This is the Aeq matrix: \n")
    print(Aeq)

    beq = np.zeros(m)
    print("This is the vector beq: \n")
    print(beq)
