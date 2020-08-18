#!/usr/bin/python3.6

import numpy as np
import scipy.sparse as sp
import gurobipy as gp
from gurobipy import GRB


if __name__ == "__main__":

    m = 2
    n = 4

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


    pre_process(A, b, Aeq, beq)

