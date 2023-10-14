# VolEsti (volume computation and sampling library)
#
# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou
#
# Contributed and/or modified by Marios Papachristou, as part of Google Summer of Code 2020 program.
#
# Licensed under GNU LGPL.3, see LICENCE file

import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sns.set_style("white")
from io import StringIO

if __name__ == '__main__':

    usage = '''
        Script to plot ODE points
        Example usage:

        ./simple_ode | python3 plot_ode.py

    '''

    argparser = argparse.ArgumentParser(usage=usage)
    argparser.add_argument('--name', type=str, default='ode_plot', help='Plot name')
    argparser.add_argument('--limits', type=str, default='auto', help='Plot limits for 2D/3D plots')
    argparser.add_argument('--step', type=float, default=1, help='Step-size (defaults to 1 if not provided)')
    args = argparser.parse_args()

    temp = StringIO(sys.stdin.read())
    data = np.loadtxt(temp)

    if len(data.shape) != 2:
        data = data.reshape(data.shape[-1], 1)

    dims = data.shape[-1]

    time = args.step * np.arange(data.shape[0])

    print('Plotting ODE components')
    fig, ax2d = plt.subplots(ncols=dims, nrows=1, squeeze=False)
    axli = ax2d.flatten()

    plt.suptitle('{}: ODE components'.format(args.name))

    for i in range(1, 1 + dims):
        axli[i-1].plot(time, data[:, i-1])
        axli[i-1].set_xlabel('$t$')
        axli[i-1].set_ylabel('$x_{}(t)$'.format(i))

    plt.savefig('{}_samples.png'.format(args.name))

    if dims in [2, 3]:
        print('Plotting trajectory (2D/3D)')
        fig = plt.figure()
        if dims == 2:
            ax = fig.add_subplot(111)
            ax.plot(data[:, 0], data[:, 1])
        if dims == 3:
            ax = fig.add_subplot(111, projection='3d')
            ax.plot(data[:, 0], data[:, 1], data[:, 2])

        ax.set_xlabel('$x_1$')
        ax.set_ylabel('$x_2$')
        if dims == 3:
            ax.set_zlabel('$x_3$')

        fig.suptitle('Trajectory'.format(args.name, dims))

        if args.limits != 'auto':
            try:
                limits = int(args.limits)
                ax.set_xlim(-limits, limits)
                ax.set_ylim(-limits, limits)
                if dims == 3:
                    ax.set_zlim(-limits, limits)
            except ValueError:
                pass

        plt.savefig('{}_trajectory.png'.format(args.name))

    plt.show()
