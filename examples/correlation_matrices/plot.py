import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# plot data
filenames = ["uniform_sampling.txt"]
for filename in filenames:
    fig = plt.figure(figsize=(4, 4))

    ax = fig.add_subplot(111, projection='3d')

    # plot points
    data = np.genfromtxt(filename, delimiter=' ')
    # ax.scatter(2,3,4) # plot the point (2,3,4) on the figure

    plt.scatter(data[:, 0], data[:, 1], data[:, 2])

# show all plots
plt.show()
