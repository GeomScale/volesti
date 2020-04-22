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

    volume = p.compute_volume(walk_len=5, epsilon=0.05, method="gaussian_annealing", seed=42)
    samples = p.generate_samples(walk_len=5, n_samples=80000, seed=42)


    print("Volume: {}".format(volume))
    print("Samples:")
    print(samples)
