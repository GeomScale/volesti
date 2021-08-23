import numpy as np
import matplotlib.pyplot as plt


A = np.array([[-1, 0],
              [0, -1],
              [1, -1],
              [-1, 1],
              [ 1, 1],
            ])
b = np.array([1,
              1,
              8,
              8,
              50])


def get_polytope_vals(x1, x2):
    X = np.vstack([x1, x2]).T
    vals = [np.max((A @ x) - b) for x in X]
    return np.array(vals)


# plot data
filenames = ["uniform_billiard_walk.txt", "uniform_accelerated_billiard_walk.txt", "gaussian_billiard_walk.txt"]
for filename in filenames:
    fig = plt.figure()

    # plot polytope
    x1_vec, x2_vec = np.meshgrid(np.linspace(-1, 30, 100), np.linspace(-1, 30, 100))
    vals = get_polytope_vals(x1_vec.ravel(), x2_vec.ravel())
    vals = vals.reshape(x1_vec.shape)
    plt.contourf(x1_vec, x2_vec, vals, levels=[-100, 0], cmap="Greys_r", alpha=0.5)

    # plot points
    data = np.genfromtxt(filename, delimiter=' ')
    plt.plot(data[:, 0], data[:, 1], 'x', c='g')
    plt.xlim(-2, 30)
    plt.ylim(-2, 30)
    plt.title(filename[:-4])

# show all plots
plt.show()