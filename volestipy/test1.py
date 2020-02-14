import numpy as np
from volestipy import HPolytope


if __name__ == "__main__":
    dim = 2
    A = np.zeros((2*dim, dim), dtype=np.float)
    A[0:dim] = np.eye(dim)
    A[dim:] -=  np.eye(dim,dim, dtype=np.float)
    b = np.ones(2*dim, dtype=np.float)


    p = HPolytope(A,b)
    print("The polytope:")
    print("Dimensions: {}".format(p.dimensions))
    print("A matrix:")
    print(p.A)
    print("b vector:")
    print(p.b)

    volume_SoB = p.compute_volume(walk_len=5, epsilon=0.05, method="sequence_of_balls", seed=42)
    volume_GA = p.compute_volume(walk_len=5, epsilon=0.05, method="gaussian_annealing", seed=42)
    samples = p.generate_samples(walk_len=5, n_samples=80000, seed=42)


    print("Volume (sequence of balls): {}".format(volume_SoB))
    print("Volume (gaussian annealing): {}".format(volume_GA))
    print("Samples:")
    print(samples)
