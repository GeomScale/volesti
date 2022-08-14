import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# plot data
filenames = ["BallWalk_matrices.txt", "RDHRWalk_matrices.txt",\
    "BilliardWalk_matrices.txt", "AcceleratedBilliardWalk_matrices.txt",\
    "BallWalk_matrices_MT.txt", "RDHRWalk_matrices_MT.txt",\
    "BilliardWalk_matrices_MT.txt", "AcceleratedBilliardWalk_matrices_MT.txt"]
for filename in filenames:
    fig = plt.figure(figsize=(4,4))
    ax = plt.axes(projection='3d')
    data = np.genfromtxt(filename, delimiter=' ')
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], s = 1)
    plt.title(filename)

# show all plots
plt.show()
