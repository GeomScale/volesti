import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# plot data
filenames = ["correlation_matrices.txt"]
for filename in filenames:
    fig = plt.figure(figsize=(4,4))
    ax = plt.axes(projection='3d')
    data = np.genfromtxt(filename, delimiter=' ')
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], s = 1)

# show all plots
plt.show()
